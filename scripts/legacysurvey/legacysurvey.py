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
from datasets import Features, Value, Array2D, Sequence, Image
from datasets.data_files import DataFilesPatternsDict
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
Image dataset from Legacy Survey DR10
"""

# TODO: Add a link to an official homepage for the dataset here
_HOMEPAGE = ""

# TODO: Add the licence for the dataset here if you can find it
_LICENSE = ""

_VERSION = "0.0.1"

_FLOAT_FEATURES = [
    "EBV",
    "FLUX_G",
    "FLUX_R",
    "FLUX_I",
    "FLUX_Z",
    "FLUX_W1",
    "FLUX_W2",
    "FLUX_W3",
    "FLUX_W4",
]

CATALOG_FEATURES = [
    "FLUX_G",
    "FLUX_R",
    "FLUX_I",
    "FLUX_Z",
    "FLUX_IVAR_G",
    "FLUX_IVAR_R",
    "FLUX_IVAR_I",
    "FLUX_IVAR_Z",
    "TYPE",
    "RA",
    "DEC",
]


class DECaLS(datasets.GeneratorBasedBuilder):
    """TODO: Short description of my dataset."""

    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="LegacySurvey",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["*/healpix=*/*.hdf5"]}
            ),
            description="LegacySurvey DR10 images.",
        ),        
        datasets.BuilderConfig(name="dr10_south_21", 
                                version=VERSION, 
                                data_files=DataFilesPatternsDict.from_patterns({'train': ['dr10_south_21/healpix=*/*.hdf5']}),
                                description="DR10 images from the southern sky, down to zmag 21"),
    ]

    DEFAULT_CONFIG_NAME = "LegacySurvey"

    _pixel_scale = 0.262
    _image_size = 160
    _bands = ['DES-G', 'DES-R', 'DES-I', 'DES-Z']

    @classmethod
    def _info(self):
        """ Defines the features available in this dataset.
        """
        # Starting with all features common to image datasets
        features = {
            "image": Sequence(
                feature={
                    "band": Value("string"),
                    "array": Array2D(
                        shape=(self._image_size, self._image_size), dtype="float32"
                    ),
                    "bit_mask": Array2D(
                        shape=(self._image_size, self._image_size), dtype="bool"
                    ),
                    "ivar": Array2D(
                        shape=(self._image_size, self._image_size), dtype="float32"
                    ),
                    "psf_fwhm": Value("float32"),
                    "scale": Value("float32"),
                }
            ),
            "model_image": Image(),
            "object_mask": Image(),
            "catalog": Sequence(
                feature={f: Value("float32") for f in CATALOG_FEATURES}
            ),
        }
        # Adding all values from the catalog
        for f in _FLOAT_FEATURES:
            features[f] = Value("float32")

        features["object_id"] = Value("string")

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
            raise ValueError(f"At least one data file must be specified, but got data_files={self.config.data_files}")
        splits = []
        for split_name, files in self.config.data_files.items():
            if isinstance(files, str):
                files = [files]
            splits.append(datasets.SplitGenerator(name=split_name, gen_kwargs={"files": files})) 
        return splits

    def _generate_examples(self, files, object_ids=None):
        """ Yields examples as (key, example) tuples.
        """
        for j, file in enumerate(files):
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                else:
                    keys = data["object_id"]
                
                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["object_id"])
                sorted_ids = data["object_id"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog 
                    i = sort_index[np.searchsorted(sorted_ids, k)]
                    # Parse image data
                    example = {
                        "image": [
                            {
                                "band": data["image_band"][i][j].decode("utf-8"),
                                "array": data["image_array"][i][j],
                                "mask": data["bit_mask"][i],
                                "ivar": data["image_ivar"][i][j],
                                "psf_fwhm": data["image_psf_fwhm"][i][j],
                                "scale": data["image_scale"][i][j],
                            }
                            for j, _ in enumerate(self._bands)
                        ],
                        "model_image": data["image_model"][i],
                        "object_mask": data["object_mask"][i],
                        "catalog": {
                            key: data[f"catalog_{key}"][i] for key in CATALOG_FEATURES
                        },
                    }
                    # Add all other requested features
                    for f in _FLOAT_FEATURES:
                        example[f] = data[f][i].astype('float32')
                    
                    # Add object type
                    example['TYPE'] = data['TYPE'][i].decode('utf-8')

                    # Add object_id
                    example["object_id"] = str(data["object_id"][i])

                    yield str(data['object_id'][i]), example