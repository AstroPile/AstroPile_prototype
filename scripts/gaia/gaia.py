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

_DESCRIPTION = """\
Spectral (BP/RP), photometric, and astrometric dataset based on Gaia DR3.
"""

_HOMEPAGE = ""

_LICENSE = ""

_VERSION = "0.0.1"

_SPECTRUM_FEATURES = [
    "coeff",
    "coeff_error",
]

_PHOTOMETRY_FEATURES = [
    "phot_g_mean_mag",
    "phot_g_mean_flux",
    "phot_g_mean_flux_error",
    "phot_bp_mean_mag",
    "phot_bp_mean_flux",
    "phot_bp_mean_flux_error",
    "phot_rp_mean_mag",
    "phot_rp_mean_flux",
    "phot_rp_mean_flux_error",
    "phot_bp_rp_excess_factor",
    "bp_rp",
    "bp_g",
    "g_rp",
]

_ASTROMETRY_FEATURES = [
    "ra",
    "ra_error",
    "dec",
    "dec_error",
    "parallax",
    "parallax_error",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "ra_dec_corr",
    "ra_parallax_corr",
    "ra_pmra_corr",
    "ra_pmdec_corr",
    "dec_parallax_corr",
    "dec_pmra_corr",
    "dec_pmdec_corr",
    "parallax_pmra_corr",
    "parallax_pmdec_corr",
    "pmra_pmdec_corr",
    "radial_velocity",
    "radial_velocity_error",
]

_GSPPHOT_FEATURES = [
    "ag_gspphot",
    "ag_gspphot_lower",
    "ag_gspphot_upper",
    "azero_gspphot",
    "azero_gspphot_lower",
    "azero_gspphot_upper",
    "distance_gspphot",
    "distance_gspphot_lower",
    "distance_gspphot_upper",
    "ebpminrp_gspphot",
    "ebpminrp_gspphot_lower",
    "ebpminrp_gspphot_upper",
    "logg_gspphot",
    "logg_gspphot_lower",
    "logg_gspphot_upper",
    "mh_gspphot",
    "mh_gspphot_lower",
    "mh_gspphot_upper",
    "teff_gspphot",
    "teff_gspphot_lower",
    "teff_gspphot_upper",
]


class Gaia(datasets.GeneratorBasedBuilder):
    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="gaia_dr3",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["gaia/healpix=*/*.hdf5"]}
            ),
            description="Gaia source table and XP coefficients.",
        ),
    ]

    DEFAULT_CONFIG_NAME = "gaia_dr3"

    @classmethod
    def _info(self):
        """Defines the features available in this dataset."""
        # Starting with all features common to image datasets
        features = {
            "spectrum": Sequence(
                {f: Value(dtype="float32") for f in _SPECTRUM_FEATURES}
            ),
            "photometry": {f: Value(dtype="float32") for f in _PHOTOMETRY_FEATURES},
            "astrometry": {f: Value(dtype="float32") for f in _ASTROMETRY_FEATURES},
            "gspphot": {f: Value(dtype="float32") for f in _GSPPHOT_FEATURES},
            "object_id": Value(dtype="int64"),
            "healpix": Value(dtype="int64"),
            "ra": Value(dtype="float32"),
            "dec": Value(dtype="float32"),
        }

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
                    keys = data["source_id"][:]

                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["source_id"][:])
                sorted_ids = data["source_id"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog
                    i = sort_index[np.searchsorted(sorted_ids, k)]

                    s_id = data["source_id"][i]

                    example = {
                        "spectrum": {f: data[f][i] for f in _SPECTRUM_FEATURES},
                        "photometry": {f: data[f][i] for f in _PHOTOMETRY_FEATURES},
                        "astrometry": {f: data[f][i] for f in _ASTROMETRY_FEATURES},
                        "gspphot": {f: data[f][i] for f in _GSPPHOT_FEATURES},
                        "object_id": s_id,
                        "healpix": data["healpix"][i],
                        "ra": data["ra"][i],
                        "dec": data["dec"][i],
                    }

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
