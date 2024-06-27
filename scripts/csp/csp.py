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
import os

_CITATION = """\
@ARTICLE{2017AJ....154..211K,
    author = {{Krisciunas}, Kevin and {Contreras}, Carlos and {Burns}, Christopher R. and {Phillips}, M.~M. and {Stritzinger}, Maximilian D. and {Morrell}, Nidia and {Hamuy}, Mario and {Anais}, Jorge and {Boldt}, Luis and {Busta}, Luis and {Campillay}, Abdo and {Castell{\'o}n}, Sergio and {Folatelli}, Gast{\'o}n and {Freedman}, Wendy L. and {Gonz{\'a}lez}, Consuelo and {Hsiao}, Eric Y. and {Krzeminski}, Wojtek and {Persson}, Sven Eric and {Roth}, Miguel and {Salgado}, Francisco and {Ser{\'o}n}, Jacqueline and {Suntzeff}, Nicholas B. and {Torres}, Sim{\'o}n and {Filippenko}, Alexei V. and {Li}, Weidong and {Madore}, Barry F. and {DePoy}, D.~L. and {Marshall}, Jennifer L. and {Rheault}, Jean-Philippe and {Villanueva}, Steven},
    title = "{The Carnegie Supernova Project. I. Third Photometry Data Release of Low-redshift Type Ia Supernovae and Other White Dwarf Explosions}",
    journal = {\aj},
    keywords = {instrumentation: photometers, supernovae: general, surveys, techniques: photometric, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
    year = 2017,
    month = nov,
    volume = {154},
    number = {5},
    eid = {211},
    pages = {211},
    doi = {10.3847/1538-3881/aa8df0},
    archivePrefix = {arXiv},
    eprint = {1709.05146},
    primaryClass = {astro-ph.IM},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2017AJ....154..211K},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

"""

_DESCRIPTION = """\
Time-series dataset from the Carnegie Supernova Project I Data Release 3 (CSP-I DR3).
"""

_HOMEPAGE = "https://csp.obs.carnegiescience.edu/"

_LICENSE = "CC BY 4.0"

_VERSION = "0.0.1"

_STR_FEATURES = [
    "object_id",
    "spec_class",
]

_FLOAT_FEATURES = [
    "ra",
    "dec",
    "redshift",
]


class CSPIDR3(datasets.GeneratorBasedBuilder):
    """"""

    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="csp_dr3",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns({"train": ["./*/*.hdf5"]}), # This seems fairly inflexible. Probably a massive failure point.
            description="Light curves from CSP-I DR3",
        ),
    ]

    DEFAULT_CONFIG_NAME = "csp_dr3"

    @classmethod
    def _info(self):
        """Defines the features available in this dataset."""
        # Starting with all features common to light curve datasets
        features = {
            'lightcurve': Sequence(feature={
                "band": Value("string"),
                "mag": Value("float32"),
                "mag_err": Value("float32"),
                "time": Value("float32"),
            }),
        }

        # Adding all values from the catalog
        for f in _FLOAT_FEATURES:
            features[f] = Value("float32")
        for f in _STR_FEATURES:
            features[f] = Value("string")

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
        for file_number, file in enumerate(itertools.chain.from_iterable(files)):
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[file_number]
                else:
                    keys = [data["object_id"][()]]

                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["object_id"][()])  # Accessing the scalar index
                sorted_ids = [data["object_id"][()]]  # Ensure this is a list of one element

                for k in keys:
                    # Extract the indices of requested ids in the catalog
                    i = sort_index[np.searchsorted(sorted_ids, k)]
                    # Parse data
                    idxs = np.arange(0, data["mag"].shape[0])
                    band_idxs = idxs.repeat(data["mag"].shape[-1]).reshape(
                        data["bands"].shape[0], -1
                    )
                    bands = [bstr.decode('utf-8') for bstr in data["bands"][()]]
                    example = {
                        'lightcurve': {
                            "band": np.asarray([bands[band_number] for band_number in band_idxs.flatten().astype("int32")]).astype("str"),
                            "time": np.asarray(data["time"]).flatten().astype("float32"),
                            "mag": np.asarray(data["mag"]).flatten().astype("float32"),
                            "mag_err": np.asarray(data["mag_err"]).flatten().astype("float32"),
                        }
                    }
                    # Add remaining features
                    for f in _FLOAT_FEATURES:
                        example[f] = np.asarray(data[f]).astype("float32")
                    for f in _STR_FEATURES:
                        # Add band names shared across dataset to each sample.
                        # I can't see a better way to do this.
                        if f == "bands":
                            example[f] = np.asarray(data[f]).astype("str")
                        else:
                            example[f] = data[f][()].astype("str")

                    yield str(data["object_id"][()]), example
