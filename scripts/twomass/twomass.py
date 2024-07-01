# Copyright 2020 The HuggingFace Datasets Authors and the current dataset script contributor.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import itertools

import numpy as np

import datasets
from datasets import Features, Value
from datasets.data_files import DataFilesPatternsDict
import h5py

_CITATION = """\
@InProceedings{huggingface:dataset,
title = {A great new dataset},
author={huggingface, Inc.
},
year={2020}
}
"""

_DESCRIPTION = """\
2MASS near-IR photometric dataset.
"""

_HOMEPAGE = ""

_LICENSE = ""

_VERSION = "0.0.1"

_mapping = dict(
    ra="float64",
    decl="float64",
    err_maj="float32",
    err_min="float32",
    err_ang="int64",
    designation="string",
    j_m="float32",
    j_cmsig="float32",
    j_msigcom="float32",
    j_snr="float32",
    h_m="float32",
    h_cmsig="float32",
    h_msigcom="float32",
    h_snr="float32",
    k_m="float32",
    k_cmsig="float32",
    k_msigcom="float32",
    k_snr="float32",
    ph_qual="string",
    rd_flg="string",
    bl_flg="string",
    cc_flg="string",
    ndet="string",
    prox="float32",
    pxpa="int64",
    pxcntr="int64",
    gal_contam="int64",
    mp_flg="int64",
    pts_key="int64",
    hemis="string",
    date="string",
    scan="int64",
    glon="float32",
    glat="float32",
    x_scan="float32",
    jdate="float64",
    j_psfchi="float32",
    h_psfchi="float32",
    k_psfchi="float32",
    j_m_stdap="float32",
    j_msig_stdap="float32",
    h_m_stdap="float32",
    h_msig_stdap="float32",
    k_m_stdap="float32",
    k_msig_stdap="float32",
    dist_edge_ns="int64",
    dist_edge_ew="int64",
    dist_edge_flg="string",
    dup_src="int64",
    use_src="int64",
    a="string",
    dist_opt="float32",
    phi_opt="int64",
    b_m_opt="float32",
    vr_m_opt="float32",
    nopt_mchs="int64",
    ext_key="int64",
    scan_key="int64",
    coadd_key="int64",
    coadd="int64",
)


class TwoMASS(datasets.GeneratorBasedBuilder):
    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="psc",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["2mass/psc/healpix=*/*.hdf5"]}
            ),
            description="2MASS point source catalog.",
        ),
    ]

    DEFAULT_CONFIG_NAME = "psc"

    @classmethod
    def _info(self):
        """Defines the features available in this dataset."""

        features = {k: Value(dtype=v, id=None) for k, v in _mapping.items()}
        features["object_id"] = Value(dtype=_mapping["pts_key"], id=None)

        return datasets.DatasetInfo(
            # This is the description that will appear on the datasets page.
            description=_DESCRIPTION,
            # This defines the different columns of the dataset and their types
            features=Features(features),
            # Homepage of the dataset for documentation
            homepage=_HOMEPAGE,
            # License for the dataset if available
            license=_LICENSE,
            # Citation for the dataset
            citation=_CITATION,
        )

    def _generate_examples(self, files, object_ids=None):
        """Yields examples as (key, example) tuples."""
        id_key = "pts_key"
        for j, file in enumerate(itertools.chain.from_iterable(files)):
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                else:
                    keys = data[id_key][:]

                # Preparing an index for fast searching through the catalog
                id_data = data[id_key][:]
                sort_index = np.argsort(id_data)
                sorted_ids = id_data[sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog
                    i = sort_index[np.searchsorted(sorted_ids, k)]

                    s_id = id_data[i]

                    example = {k: data[k][i] for k in _mapping}
                    example["object_id"] = s_id

                    for k, v in example.items():
                        if isinstance(v, bytes):
                            example[k] = v.decode("utf-8")
                        if isinstance(v, float) and np.isnan(v):
                            example[k] = None

                    yield int(s_id), example

    def _split_generators(self, dl_manager):
        """We handle string, list and dicts in datafiles"""
        if not self.config.data_files:
            raise ValueError(
                f"At least one data file must be specified, but got data_files={self.config.data_files}"
            )
        data_files = dl_manager.download_and_extract(self.config.data_files)
        if isinstance(data_files, (str, list, tuple)):
            files = data_files
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            return [
                datasets.SplitGenerator(
                    name=datasets.Split.TRAIN, gen_kwargs={"files": files}
                )
            ]
        splits = []
        for split_name, files in data_files.items():
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            splits.append(
                datasets.SplitGenerator(name=split_name, gen_kwargs={"files": files})
            )
        return splits
