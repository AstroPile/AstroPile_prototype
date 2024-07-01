from typing import List, Optional, Tuple, Union

from transformers import InformerPreTrainedModel, InformerModel, InformerForPrediction, InformerConfig
from transformers.models.informer.modeling_informer import InformerConvLayer, InformerEncoder, InformerEncoderLayer, InformerDecoder, _prepare_4d_attention_mask, weighted_average, nll
from transformers.modeling_outputs import SequenceClassifierOutput, Seq2SeqTSModelOutput, BaseModelOutput, Seq2SeqTSPredictionOutput, BaseModelOutputWithPastAndCrossAttentions, MaskedLMOutput
from transformers.time_series_utils import StudentTOutput, NormalOutput, NegativeBinomialOutput

import torch
import torch.nn as nn
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, MSELoss, functional as F

from connect_later.dataset_preprocess_raw import create_network_inputs

import pdb
import numpy as np
import random

class MultiDimFourierPE(nn.Module):
    """
    Implements Li et al. 2021 learnable multi-dimensional fourier positional encoding
    """

    def __init__(self, fourier_dim: int, hidden_dim: int, model_dim: int) -> None:
        super().__init__()
        self.fourier_dim = fourier_dim # e.g. 384
        self.hidden_dim = hidden_dim # e.g. 32
        self.embed_dim = model_dim # e.g. 768

        # initialize layers here to avoid device incompatibilities
        self.fourier_to_hidden = nn.Linear(self.fourier_dim, self.hidden_dim)
        self.hidden_to_embed = nn.Linear(self.hidden_dim, self.embed_dim)

    def forward(self, pos: torch.Tensor) -> torch.Tensor:
        """
        Args:
            pos: bsz x seq_len x num_dimensions
            - num_dimensions should be 2: [:,0] = time, [:,1] = central wavelength of band
            - both already scaled to be between 0 and 1

        Returns:
            pos_encoding: bsz x seq_len x d_model
        """
        bsz, seq_len, num_dimensions = pos.size()
        out = torch.zeros(bsz, seq_len, self.fourier_dim, device=pos.device, dtype=pos.dtype)
        sentinel = self.fourier_dim // 2 if self.fourier_dim % 2 == 0 else (self.fourier_dim // 2) + 1

        fourier_coeff = nn.Linear(num_dimensions, sentinel, bias=False, device=pos.device) # map num_dimensions to half of the fourier feature dimension
        out[:, :, 0:sentinel] = torch.sin(fourier_coeff(pos))
        out[:, :, sentinel:] = torch.cos(fourier_coeff(pos))
        out *= (1 / np.sqrt(self.fourier_dim)) # scale by sqrt of fourier dimension

        out = self.fourier_to_hidden(out) # map to hidden dimension
        out = nn.GELU()(out)
        out = self.hidden_to_embed(out) # map to embedding dimension

        return out

class InformerEncoderFourierPE(InformerEncoder):
    """
    Informer encoder consisting of *config.encoder_layers* self attention layers with distillation layers. Each
    attention layer is an [`InformerEncoderLayer`]. This implementation includes a custom positional encoding, CustomPE.

    Args:
        config: InformerConfig
    """

    def __init__(self, config: InformerConfig):
        super().__init__(config)
        print("Using Fourier PE")
        fourier_dim = config.fourier_dim if hasattr(config, "fourier_dim") else 384
        hidden_dim = config.PE_hidden_dim if hasattr(config, "PE_hidden_dim") else 32
        self.embed_positions = MultiDimFourierPE(fourier_dim, hidden_dim, config.d_model)

        # Initialize weights and apply final processing
        self.post_init()

    # copied from InformerEncoder
    def forward(
        self,
        attention_mask: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.FloatTensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple, BaseModelOutput]:
        r"""
        Args:
            attention_mask (`torch.Tensor` of shape `(batch_size, sequence_length)`, *optional*):
                Mask to avoid performing attention on padding token indices. Mask values selected in `[0, 1]`:

                - 1 for tokens that are **not masked**,
                - 0 for tokens that are **masked**.

                [What are attention masks?](../glossary#attention-mask)
            head_mask (`torch.Tensor` of shape `(encoder_layers, encoder_attention_heads)`, *optional*):
                Mask to nullify selected heads of the attention modules. Mask values selected in `[0, 1]`:

                - 1 indicates the head is **not masked**,
                - 0 indicates the head is **masked**.

            inputs_embeds (`torch.FloatTensor` of shape `(batch_size, sequence_length, hidden_size)`, *optional*):
                Optionally, instead of passing `input_ids` you can choose to directly pass an embedded representation.
                This is useful if you want more control over how to convert `input_ids` indices into associated vectors
                than the model's internal embedding lookup matrix.
            output_attentions (`bool`, *optional*):
                Whether or not to return the attentions tensors of all attention layers. See `attentions` under
                returned tensors for more detail.
            output_hidden_states (`bool`, *optional*):
                Whether or not to return the hidden states of all layers. See `hidden_states` under returned tensors
                for more detail.
            return_dict (`bool`, *optional*):
                Whether or not to return a [`~utils.ModelOutput`] instead of a plain tuple.
        """
        output_attentions = output_attentions if output_attentions is not None else self.config.output_attentions
        output_hidden_states = (
            output_hidden_states if output_hidden_states is not None else self.config.output_hidden_states
        )
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict

        # inputs_embeds[:,:,0]:
            # flux     XXXX
            # flux_err XXXX
            # mean     XXXX from scaling, not sure why they have 2 rows each
            # mean     XXXX
            # std      XXXX
            # std      XXXX
            # time     XXXX
            # wv       XXXX
        hidden_states = self.value_embedding(inputs_embeds)
        embed_pos = self.embed_positions(inputs_embeds[:,:,-2:]) # last two rows of features dimension are time and central wavelength

        hidden_states = self.layernorm_embedding(hidden_states + embed_pos)
        hidden_states = nn.functional.dropout(hidden_states, p=self.dropout, training=self.training)

        # expand attention_mask
        if attention_mask is not None:
            # [bsz, seq_len] -> [bsz, 1, tgt_seq_len, src_seq_len]
            attention_mask = _prepare_4d_attention_mask(attention_mask, inputs_embeds.dtype)

        encoder_states = () if output_hidden_states else None
        all_attentions = () if output_attentions else None

        # check if head_mask has a correct number of layers specified if desired
        if head_mask is not None:
            if head_mask.size()[0] != (len(self.layers)):
                raise ValueError(
                    f"The head_mask should be specified for {len(self.layers)} layers, but it is for"
                    f" {head_mask.size()[0]}."
                )

        for idx, (encoder_layer, conv_layer) in enumerate(zip(self.layers, self.conv_layers)):
            if output_hidden_states:
                encoder_states = encoder_states + (hidden_states,)
            # add LayerDrop (see https://arxiv.org/abs/1909.11556 for description)
            dropout_probability = random.uniform(0, 1)
            if self.training and (dropout_probability < self.layerdrop):  # skip the layer
                layer_outputs = (None, None)
            else:
                if self.gradient_checkpointing and self.training:

                    def create_custom_forward(module):
                        def custom_forward(*inputs):
                            return module(*inputs, output_attentions)

                        return custom_forward

                    layer_outputs = torch.utils.checkpoint.checkpoint(
                        create_custom_forward(encoder_layer),
                        hidden_states,
                        attention_mask,
                        (head_mask[idx] if head_mask is not None else None),
                    )
                    if conv_layer is not None:
                        output = torch.utils.checkpoint.checkpoint(conv_layer, layer_outputs[0])
                        layer_outputs = (output,) + layer_outputs[1:]
                else:
                    layer_outputs = encoder_layer(
                        hidden_states,
                        attention_mask,
                        layer_head_mask=(head_mask[idx] if head_mask is not None else None),
                        output_attentions=output_attentions,
                    )
                    if conv_layer is not None:
                        output = conv_layer(layer_outputs[0])
                        layer_outputs = (output,) + layer_outputs[1:]

                hidden_states = layer_outputs[0]

            if output_attentions:
                all_attentions = all_attentions + (layer_outputs[1],)

        if output_hidden_states:
            encoder_states = encoder_states + (hidden_states,)

        if not return_dict:
            return tuple(v for v in [hidden_states, encoder_states, all_attentions] if v is not None)
        return BaseModelOutput(
            last_hidden_state=hidden_states, hidden_states=encoder_states, attentions=all_attentions
        )

class InformerForSequenceClassification(InformerPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)
        self.num_labels = config.num_labels
        print(f"num labels: {self.num_labels}")
        self.config = config

        # if config.fourier_pe and config.mask:
        config.distil = False
        self.encoder = InformerEncoderFourierPE(config)
        self.dropout = nn.Dropout(0.3)
        self.classifier = nn.Linear(config.hidden_size, config.num_labels)

        # Initialize weights and apply final processing
        self.post_init()

    def forward(
        self,
        past_values: torch.Tensor,
        past_time_features: torch.Tensor,
        past_observed_mask: torch.Tensor,
        labels: Optional[torch.Tensor] = None,
        weights: Optional[torch.Tensor] = None,
        static_categorical_features: Optional[torch.Tensor] = None,
        static_real_features: Optional[torch.Tensor] = None,
        future_values: Optional[torch.Tensor] = None,
        future_time_features: Optional[torch.Tensor] = None,
        future_observed_mask: Optional[torch.Tensor] = None,
        decoder_attention_mask: Optional[torch.LongTensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        decoder_head_mask: Optional[torch.Tensor] = None,
        cross_attn_head_mask: Optional[torch.Tensor] = None,
        encoder_outputs: Optional[List[torch.FloatTensor]] = None,
        past_key_values: Optional[List[torch.FloatTensor]] = None,
        output_hidden_states: Optional[bool] = None,
        output_attentions: Optional[bool] = None,
        use_cache: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[SequenceClassifierOutput, Tuple]:

        return_dict = return_dict if return_dict is not None else self.config.use_return_dict

        if self.config.mask:
            #TODO: should encapsulate this preprocessing + encoder into its own class
            transformer_inputs, loc, scale, static_feat = create_network_inputs(
                config=self.config,
                past_values=past_values,
                past_time_features=past_time_features,
                past_observed_mask=past_observed_mask,
                static_categorical_features=static_categorical_features,
                static_real_features=static_real_features,
            )

            outputs = self.encoder(
                inputs_embeds=transformer_inputs,
                head_mask=head_mask,
                output_attentions=output_attentions,
                output_hidden_states=output_hidden_states,
                return_dict=return_dict,
            )

            decoder_output = outputs.last_hidden_state
            pooled_output = torch.mean(decoder_output, dim=1, keepdim=True) # average over time dimension

        else:
            outputs = self.model(
                past_values=past_values,
                past_time_features=past_time_features,
                past_observed_mask=past_observed_mask,
                future_values=future_values,
                future_time_features=future_time_features,
                static_categorical_features=static_categorical_features,
                static_real_features=static_real_features,
                decoder_attention_mask=decoder_attention_mask,
                head_mask=head_mask,
                decoder_head_mask=decoder_head_mask,
                cross_attn_head_mask=cross_attn_head_mask,
                encoder_outputs=encoder_outputs,
                past_key_values=past_key_values,
                output_hidden_states=output_hidden_states,
                output_attentions=output_attentions,
                use_cache=use_cache,
                return_dict=return_dict,
            )

            decoder_output = outputs.last_hidden_state

            pooled_output = torch.mean(decoder_output, dim=1, keepdim=True) # average over time dimension
            # for _ in range(self.num_conv_layers):
            #     decoder_output = self.conv_layer(decoder_output)
            # pooled_output = self.pooler_activation(self.pooler(decoder_output))

        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)

        loss = None
        if labels is not None:
            if self.config.num_labels == 1 and self.config.regression:
                loss_fn = MSELoss()
            elif self.config.num_labels == 1:
                labels = labels.float()
                loss_fn = BCEWithLogitsLoss()
            else:
                loss_fn = CrossEntropyLoss(weight=weights) if weights is not None else CrossEntropyLoss()
            loss = loss_fn(torch.squeeze(logits, 1), labels) # avoid squeezing batch dimension if batch has 1 object

        if not return_dict:
            output = (logits,) + outputs[2:]
            return ((loss,) + output) if loss is not None else output

        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            # hidden_states=outputs.encoder_hidden_states,
            hidden_states=outputs.last_hidden_state,
            # attentions=outputs.encoder_attentions,
        )
