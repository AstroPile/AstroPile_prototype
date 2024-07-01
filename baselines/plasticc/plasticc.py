import argparse
from pathlib import Path
import yaml

import torch
from transformers import InformerConfig
from datasets import load_dataset, load_metric

from models import InformerForSequenceClassification
from preprocess import create_test_dataloader_raw

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', type=str)
    parser.add_argument('--pretrained_model_path', type=str)
    return parser.parse_args()

def get_config(config):

    model_config = InformerConfig(
        input_size=2,
        prediction_length=0,
        context_length=300,
        lags_sequence=[0],
        num_time_features=2, #wavelength + time
        num_static_real_features=0 ,

        # informer params:
        dropout=config['dropout_rate'],
        encoder_layers=config['num_encoder_layers'],
        decoder_layers=config['num_decoder_layers'],
        d_model=config['d_model'],
        scaling=None,
        has_labels=False,
        mask=True,
        mask_probability=args.mask_probability,
    )

    return model_config

def main():
    args = parse_args()
    config = yaml.safe_load("plasticc_config.yml")
    model_config = get_config(config)

    model = InformerForSequenceClassification.from_pretrained(
           args.pretrained_model_path,
           config=model_config,
           ignore_mismatched_sizes=True
    )

    dataset = load_dataset(args.data_path)
    val_dataloader = create_test_dataloader_raw(
        dataset=dataset,
        batch_size=config["batch_size"],
        compute_loss=True,# no longer optional for encoder-decoder latent space
        has_labels=True,
    )

    model.eval()
    metric = load_metric("accuracy")
    for idx, batch in enumerate(val_dataloader):
        with torch.no_grad():
            outputs = model(**batch)
        metric.add_batch(predictions=outputs.logits, references=batch["labels"])
    print(f"Accuracy: {metric.compute()}")
