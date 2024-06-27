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
AllWISE mid-IR photometric dataset.
"""

_HOMEPAGE = ""

_LICENSE = ""

_VERSION = "0.0.1"

_mapping = dict(
    designation="string",
    ra="float32",
    dec="float32",
    sigra="float32",
    sigdec="float32",
    sigradec="float32",
    glon="float32",
    glat="float32",
    elon="float32",
    elat="float32",
    wx="float32",
    wy="float32",
    cntr="int64",
    source_id="string",
    coadd_id="string",
    src="int64",
    w1mpro="float32",
    w1sigmpro="float32",
    w1snr="float32",
    w1rchi2="float32",
    w2mpro="float32",
    w2sigmpro="float32",
    w2snr="float32",
    w2rchi2="float32",
    w3mpro="float32",
    w3sigmpro="float32",
    w3snr="float32",
    w3rchi2="float32",
    w4mpro="float32",
    w4sigmpro="float32",
    w4snr="float32",
    w4rchi2="float32",
    rchi2="float32",
    nb="int64",
    na="int64",
    w1sat="float32",
    w2sat="float32",
    w3sat="float32",
    w4sat="float32",
    satnum="string",
    ra_pm="float32",
    dec_pm="float32",
    sigra_pm="float32",
    sigdec_pm="float32",
    sigradec_pm="float32",
    pmra="int64",
    sigpmra="int64",
    pmdec="int64",
    sigpmdec="int64",
    w1rchi2_pm="float32",
    w2rchi2_pm="float32",
    w3rchi2_pm="float32",
    w4rchi2_pm="float32",
    rchi2_pm="float32",
    pmcode="string",
    cc_flags="string",
    rel="string",
    ext_flg="int64",
    var_flg="string",
    ph_qual="string",
    det_bit="int64",
    moon_lev="string",
    w1nm="int64",
    w1m="int64",
    w2nm="int64",
    w2m="int64",
    w3nm="int64",
    w3m="int64",
    w4nm="int64",
    w4m="int64",
    w1cov="float32",
    w2cov="float32",
    w3cov="float32",
    w4cov="float32",
    w1cc_map="int64",
    w1cc_map_str="string",
    w2cc_map="int64",
    w2cc_map_str="string",
    w3cc_map="int64",
    w3cc_map_str="string",
    w4cc_map="int64",
    w4cc_map_str="string",
    best_use_cntr="int64",
    ngrp="int64",
    w1flux="float32",
    w1sigflux="float32",
    w1sky="float32",
    w1sigsk="float32",
    w1conf="float32",
    w2flux="float32",
    w2sigflux="float32",
    w2sky="float32",
    w2sigsk="float32",
    w2conf="float32",
    w3flux="float32",
    w3sigflux="float32",
    w3sky="float32",
    w3sigsk="float32",
    w3conf="float32",
    w4flux="float32",
    w4sigflux="float32",
    w4sky="float32",
    w4sigsk="float32",
    w4conf="float32",
    w1mag="float32",
    w1sigm="float32",
    w1flg="int64",
    w1mcor="float32",
    w2mag="float32",
    w2sigm="float32",
    w2flg="int64",
    w2mcor="float32",
    w3mag="float32",
    w3sigm="float32",
    w3flg="int64",
    w3mcor="float32",
    w4mag="float32",
    w4sigm="float32",
    w4flg="int64",
    w4mcor="float32",
    w1mag_1="float32",
    w1sigm_1="float32",
    w1flg_1="int64",
    w2mag_1="float32",
    w2sigm_1="float32",
    w2flg_1="int64",
    w3mag_1="float32",
    w3sigm_1="float32",
    w3flg_1="int64",
    w4mag_1="float32",
    w4sigm_1="float32",
    w4flg_1="int64",
    w1mag_2="float32",
    w1sigm_2="float32",
    w1flg_2="int64",
    w2mag_2="float32",
    w2sigm_2="float32",
    w2flg_2="int64",
    w3mag_2="float32",
    w3sigm_2="float32",
    w3flg_2="int64",
    w4mag_2="float32",
    w4sigm_2="float32",
    w4flg_2="int64",
    w1mag_3="float32",
    w1sigm_3="float32",
    w1flg_3="int64",
    w2mag_3="float32",
    w2sigm_3="float32",
    w2flg_3="int64",
    w3mag_3="float32",
    w3sigm_3="float32",
    w3flg_3="int64",
    w4mag_3="float32",
    w4sigm_3="float32",
    w4flg_3="int64",
    w1mag_4="float32",
    w1sigm_4="float32",
    w1flg_4="int64",
    w2mag_4="float32",
    w2sigm_4="float32",
    w2flg_4="int64",
    w3mag_4="float32",
    w3sigm_4="float32",
    w3flg_4="int64",
    w4mag_4="float32",
    w4sigm_4="float32",
    w4flg_4="int64",
    w1mag_5="float32",
    w1sigm_5="float32",
    w1flg_5="int64",
    w2mag_5="float32",
    w2sigm_5="float32",
    w2flg_5="int64",
    w3mag_5="float32",
    w3sigm_5="float32",
    w3flg_5="int64",
    w4mag_5="float32",
    w4sigm_5="float32",
    w4flg_5="int64",
    w1mag_6="float32",
    w1sigm_6="float32",
    w1flg_6="int64",
    w2mag_6="float32",
    w2sigm_6="float32",
    w2flg_6="int64",
    w3mag_6="float32",
    w3sigm_6="float32",
    w3flg_6="int64",
    w4mag_6="float32",
    w4sigm_6="float32",
    w4flg_6="int64",
    w1mag_7="float32",
    w1sigm_7="float32",
    w1flg_7="int64",
    w2mag_7="float32",
    w2sigm_7="float32",
    w2flg_7="int64",
    w3mag_7="float32",
    w3sigm_7="float32",
    w3flg_7="int64",
    w4mag_7="float32",
    w4sigm_7="float32",
    w4flg_7="int64",
    w1mag_8="float32",
    w1sigm_8="float32",
    w1flg_8="int64",
    w2mag_8="float32",
    w2sigm_8="float32",
    w2flg_8="int64",
    w3mag_8="float32",
    w3sigm_8="float32",
    w3flg_8="int64",
    w4mag_8="float32",
    w4sigm_8="float32",
    w4flg_8="int64",
    w1magp="float32",
    w1sigp1="float32",
    w1sigp2="float32",
    w1k="float32",
    w1ndf="int64",
    w1mlq="float32",
    w1mjdmin="float32",
    w1mjdmax="float32",
    w1mjdmean="float32",
    w2magp="float32",
    w2sigp1="float32",
    w2sigp2="float32",
    w2k="float32",
    w2ndf="int64",
    w2mlq="float32",
    w2mjdmin="float32",
    w2mjdmax="float32",
    w2mjdmean="float32",
    w3magp="float32",
    w3sigp1="float32",
    w3sigp2="float32",
    w3k="float32",
    w3ndf="int64",
    w3mlq="float32",
    w3mjdmin="float32",
    w3mjdmax="float32",
    w3mjdmean="float32",
    w4magp="float32",
    w4sigp1="float32",
    w4sigp2="float32",
    w4k="float32",
    w4ndf="int64",
    w4mlq="float32",
    w4mjdmin="float32",
    w4mjdmax="float32",
    w4mjdmean="float32",
    rho12="int64",
    rho23="int64",
    rho34="int64",
    q12="int64",
    q23="int64",
    q34="int64",
    xscprox="float32",
    w1rsemi="float32",
    w1ba="float32",
    w1pa="float32",
    w1gmag="float32",
    w1gerr="float32",
    w1gflg="int64",
    w2rsemi="float32",
    w2ba="float32",
    w2pa="float32",
    w2gmag="float32",
    w2gerr="float32",
    w2gflg="int64",
    w3rsemi="float32",
    w3ba="float32",
    w3pa="float32",
    w3gmag="float32",
    w3gerr="float32",
    w3gflg="int64",
    w4rsemi="float32",
    w4ba="float32",
    w4pa="float32",
    w4gmag="float32",
    w4gerr="float32",
    w4gflg="int64",
    tmass_key="int64",
    r_2mass="float32",
    pa_2mass="float32",
    n_2mass="int64",
    j_m_2mass="float32",
    j_msig_2mass="float32",
    h_m_2mass="float32",
    h_msig_2mass="float32",
    k_m_2mass="float32",
    k_msig_2mass="float32",
    x="float32",
    y="float32",
    z="float32",
    spt_ind="int64",
    htm20="int64",
)


class AllWISE(datasets.GeneratorBasedBuilder):
    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="allwise",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["allwise/healpix=*/*.hdf5"]}
            ),
            description="AllWISE source catalog.",
        ),
    ]

    DEFAULT_CONFIG_NAME = "allwise"

    @classmethod
    def _info(self):
        """Defines the features available in this dataset."""

        features = {k: Value(dtype=v, id=None) for k, v in _mapping.items()}
        features["object_id"] = Value(dtype=_mapping["cntr"], id=None)

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
        for j, file in enumerate(itertools.chain.from_iterable(files)):
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                else:
                    keys = data["cntr"][:]

                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["cntr"][:])
                sorted_ids = data["cntr"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog
                    i = sort_index[np.searchsorted(sorted_ids, k)]

                    s_id = data["cntr"][i]

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
