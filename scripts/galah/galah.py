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
import datasets
from datasets import Features, Value, Sequence
from datasets.data_files import DataFilesPatternsDict
import itertools
import h5py
import numpy as np

# TODO: Add BibTeX citation
# Find for instance the citation on arxiv or on the dataset repo/website
_CITATION = """\
@InProceedings{huggingface:dataset,
title = {A great new dataset},
author={huggingface, Inc.
},
year={2020}
}
"""

# TODO: Add description of the dataset here
# You can copy an official description
_DESCRIPTION = """\
Spectra dataset based on GALAH.
"""

# TODO: Add a link to an official homepage for the dataset here
_HOMEPAGE = ""

# TODO: Add the licence for the dataset here if you can find it
_LICENSE = ""

_VERSION = "0.0.1"

# Full list of features available here:
# https://data.sdss.org/datamodel/files/SPECTRO_REDUX/specObj.html
_FLOAT_FEATURES = [
    'timestamp',
    'ra',
    'dec',
    'teff',
    'e_teff',
    'logg',
    'e_logg',
    'fe_h',
    'e_fe_h',
    'fe_h_atmo',
    'vmic',
    'vbroad',
    'e_vbroad',
    'alpha_fe',
    'e_alpha_fe'
]

class GALAH(datasets.GeneratorBasedBuilder):

    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="galah_dr3",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["galah_dr3/healpix=*/*.hdf5"]}
            ),
            description="GALAH DR3",
        ),
    ]

    DEFAULT_CONFIG_NAME = "galah_dr3"

    @classmethod
    def _info(self):
        """Defines the features available in this dataset."""
        # Starting with all features common to image datasets
        features = {
            "spectrum": Sequence({
                "flux": Value(dtype="float32"),
                "ivar": Value(dtype="float32"),
                "lsf": Value(dtype="float32"),
                "lsf_sigma":  Value(dtype="float32"),
                "lambda": Value(dtype="float32"),
                "norm_flux": Value(dtype="float32"),
                "norm_ivar": Value(dtype="float32"),
                "norm_lambda": Value(dtype="float32"),
            }),
            "filter_indices": {
                "B_start": Value(dtype="int32"),
                "B_end": Value(dtype="int32"),
                "G_start": Value(dtype="int32"),
                "G_end": Value(dtype="int32"),
                "R_start": Value(dtype="int32"),
                "R_end": Value(dtype="int32"),
                "I_start": Value(dtype="int32"),
                "I_end": Value(dtype="int32")
            }
        }

        # Adding all values from the catalog
        for f in _FLOAT_FEATURES:
            features[f] = Value("float32")
        
        features["object_id"] = Value("string")
        features["flux_unit"] = Value("string")

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

    def _generate_examples(self, files, object_ids=None):
        """Yields examples as (key, example) tuples."""
        for j, file in enumerate(itertools.chain.from_iterable(files)):
            with h5py.File(file, "r+") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                else:
                    keys = data["object_id"][:]

                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["object_id"][:])
                sorted_ids = data["object_id"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog
                    i = sort_index[np.searchsorted(sorted_ids, k)]

                    # Parse spectrum data
                    example = {
                        "spectrum": {
                            "flux": data["spectrum_flux"][i],
                            "ivar": data["spectrum_flux_ivar"][i],
                            "lsf_sigma": data["spectrum_lsf_sigma"][i],
                            "lambda": data["spectrum_lambda"][i],
                            "norm_flux": data["spectrum_norm_flux"][i],
                            "norm_ivar": data["spectrum_norm_ivar"][i],
                            "norm_lambda": data["spectrum_norm_lambda"][i]
                        }
                    }
                    example["filter_indices"] = {
                        "B_start": data["spectrum_B_ind_start"][i].astype("int32"),
                        "B_end": data["spectrum_B_ind_end"][i].astype("int32"),
                        "G_start": data["spectrum_G_ind_start"][i].astype("int32"),
                        "G_end": data["spectrum_G_ind_end"][i].astype("int32"),
                        "R_start": data["spectrum_R_ind_start"][i].astype("int32"),
                        "R_end": data["spectrum_R_ind_end"][i].astype("int32"),
                        "I_start": data["spectrum_I_ind_start"][i].astype("int32"),
                        "I_end": data["spectrum_I_ind_end"][i].astype("int32")
                    }
                    # Add all other requested features
                    for f in _FLOAT_FEATURES:
                        example[f] = data[f][i].astype("float32")

                    # Add object_id
                    example["object_id"] = str(data["object_id"][i])
                    example["flux_unit"] = "erg/cm^2/A/s"

                    yield str(data["object_id"][i]), example
