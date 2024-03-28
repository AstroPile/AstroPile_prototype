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
from datasets import Features, Value, Array2D, Sequence
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
Image dataset based on HSC SSP PRD3.
"""

# TODO: Add a link to an official homepage for the dataset here
_HOMEPAGE = ""

# TODO: Add the licence for the dataset here if you can find it
_LICENSE = ""

_VERSION = "0.0.1"

_FLOAT_FEATURES = [
    'thresh',
'npix', 'tnpix', 'xmin', 'xmax', 'ymin', 'ymax', 'x', 'y', 'x2_image', 'y2_image', 'xy_image', 'errx2', 'erry2', 'errxy', 'a_image', 'b_image', 'theta_image', 'cxx_image', 'cyy_image', 'cxy_image', 'cflux', 'flux', 'cpeak', 'peak', 'xcpeak', 'ycpeak', 'xpeak', 'ypeak', 'flag', 'x_image', 'y_image', 'number', 'ra', 'dec', 'x_world', 'y_world', 'flux_iso', 'fluxerr_iso', 'area_iso', 'mag_iso', 'kron_radius', 'kron_rcirc', 'flux_auto', 'fluxerr_auto', 'bkg_auto', 'flag_auto', 'area_auto', 'flux_radius_flag', 'flux_radius_20', 'flux_radius', 'flux_radius_90', 'tot_corr', 'mag_auto', 'magerr_auto', 'flux_aper_0', 'fluxerr_aper_0', 'flag_aper_0', 'bkg_aper_0', 'mask_aper_0', 'flux_aper_1', 'fluxerr_aper_1', 'flag_aper_1', 'bkg_aper_1', 'mask_aper_1', 'flux_aper_2', 'fluxerr_aper_2', 'flag_aper_2', 'bkg_aper_2', 'mask_aper_2', 'f105w_flux_aper_0', 'f105w_fluxerr_aper_0', 'f105w_flag_aper_0', 'f105w_bkg_aper_0', 'f105w_mask_aper_0', 'f105w_flux_aper_1', 'f105w_fluxerr_aper_1', 'f105w_flag_aper_1', 'f105w_bkg_aper_1', 'f105w_mask_aper_1', 'f105w_flux_aper_2', 'f105w_fluxerr_aper_2', 'f105w_flag_aper_2', 'f105w_bkg_aper_2', 'f105w_mask_aper_2', 'f105w_tot_corr', 'f115w_flux_aper_0', 'f115w_fluxerr_aper_0', 'f115w_flag_aper_0', 'f115w_bkg_aper_0', 'f115w_mask_aper_0', 'f115w_flux_aper_1', 'f115w_fluxerr_aper_1', 'f115w_flag_aper_1', 'f115w_bkg_aper_1', 'f115w_mask_aper_1', 'f115w_flux_aper_2', 'f115w_fluxerr_aper_2', 'f115w_flag_aper_2', 'f115w_bkg_aper_2', 'f115w_mask_aper_2', 'f115w_tot_corr', 'f125w_flux_aper_0', 'f125w_fluxerr_aper_0', 'f125w_flag_aper_0', 'f125w_bkg_aper_0', 'f125w_mask_aper_0', 'f125w_flux_aper_1', 'f125w_fluxerr_aper_1', 'f125w_flag_aper_1', 'f125w_bkg_aper_1', 'f125w_mask_aper_1', 'f125w_flux_aper_2', 'f125w_fluxerr_aper_2', 'f125w_flag_aper_2', 'f125w_bkg_aper_2', 'f125w_mask_aper_2', 'f125w_tot_corr', 'f140w_flux_aper_0', 'f140w_fluxerr_aper_0', 'f140w_flag_aper_0', 'f140w_bkg_aper_0', 'f140w_mask_aper_0', 'f140w_flux_aper_1', 'f140w_fluxerr_aper_1', 'f140w_flag_aper_1', 'f140w_bkg_aper_1', 'f140w_mask_aper_1', 'f140w_flux_aper_2', 'f140w_fluxerr_aper_2', 'f140w_flag_aper_2', 'f140w_bkg_aper_2', 'f140w_mask_aper_2', 'f140w_tot_corr', 'f150w_flux_aper_0', 'f150w_fluxerr_aper_0', 'f150w_flag_aper_0', 'f150w_bkg_aper_0', 'f150w_mask_aper_0', 'f150w_flux_aper_1', 'f150w_fluxerr_aper_1', 'f150w_flag_aper_1', 'f150w_bkg_aper_1', 'f150w_mask_aper_1', 'f150w_flux_aper_2', 'f150w_fluxerr_aper_2', 'f150w_flag_aper_2', 'f150w_bkg_aper_2', 'f150w_mask_aper_2', 'f150w_tot_corr', 'f160w_flux_aper_0', 'f160w_fluxerr_aper_0', 'f160w_flag_aper_0', 'f160w_bkg_aper_0', 'f160w_mask_aper_0', 'f160w_flux_aper_1', 'f160w_fluxerr_aper_1', 'f160w_flag_aper_1', 'f160w_bkg_aper_1', 'f160w_mask_aper_1', 'f160w_flux_aper_2', 'f160w_fluxerr_aper_2', 'f160w_flag_aper_2', 'f160w_bkg_aper_2', 'f160w_mask_aper_2', 'f160w_tot_corr', 'f200w_flux_aper_0', 'f200w_fluxerr_aper_0', 'f200w_flag_aper_0', 'f200w_bkg_aper_0', 'f200w_mask_aper_0', 'f200w_flux_aper_1', 'f200w_fluxerr_aper_1', 'f200w_flag_aper_1', 'f200w_bkg_aper_1', 'f200w_mask_aper_1', 'f200w_flux_aper_2', 'f200w_fluxerr_aper_2', 'f200w_flag_aper_2', 'f200w_bkg_aper_2', 'f200w_mask_aper_2', 'f200w_tot_corr', 'f277w_flux_aper_0', 'f277w_fluxerr_aper_0', 'f277w_flag_aper_0', 'f277w_bkg_aper_0', 'f277w_mask_aper_0', 'f277w_flux_aper_1', 'f277w_fluxerr_aper_1', 'f277w_flag_aper_1', 'f277w_bkg_aper_1', 'f277w_mask_aper_1', 'f277w_flux_aper_2', 'f277w_fluxerr_aper_2', 'f277w_flag_aper_2', 'f277w_bkg_aper_2', 'f277w_mask_aper_2', 'f277w_tot_corr', 'f356w_flux_aper_0', 'f356w_fluxerr_aper_0', 'f356w_flag_aper_0', 'f356w_bkg_aper_0', 'f356w_mask_aper_0', 'f356w_flux_aper_1', 'f356w_fluxerr_aper_1', 'f356w_flag_aper_1', 'f356w_bkg_aper_1', 'f356w_mask_aper_1', 'f356w_flux_aper_2', 'f356w_fluxerr_aper_2', 'f356w_flag_aper_2', 'f356w_bkg_aper_2', 'f356w_mask_aper_2', 'f356w_tot_corr', 'f410m_flux_aper_0', 'f410m_fluxerr_aper_0', 'f410m_flag_aper_0', 'f410m_bkg_aper_0', 'f410m_mask_aper_0', 'f410m_flux_aper_1', 'f410m_fluxerr_aper_1', 'f410m_flag_aper_1', 'f410m_bkg_aper_1', 'f410m_mask_aper_1', 'f410m_flux_aper_2', 'f410m_fluxerr_aper_2', 'f410m_flag_aper_2', 'f410m_bkg_aper_2', 'f410m_mask_aper_2', 'f410m_tot_corr', 'f435w_flux_aper_0', 'f435w_fluxerr_aper_0', 'f435w_flag_aper_0', 'f435w_bkg_aper_0', 'f435w_mask_aper_0', 'f435w_flux_aper_1', 'f435w_fluxerr_aper_1', 'f435w_flag_aper_1', 'f435w_bkg_aper_1', 'f435w_mask_aper_1', 'f435w_flux_aper_2', 'f435w_fluxerr_aper_2', 'f435w_flag_aper_2', 'f435w_bkg_aper_2', 'f435w_mask_aper_2', 'f435w_tot_corr', 'f444w_flux_aper_0', 'f444w_fluxerr_aper_0', 'f444w_flag_aper_0', 'f444w_bkg_aper_0', 'f444w_mask_aper_0', 'f444w_flux_aper_1', 'f444w_fluxerr_aper_1', 'f444w_flag_aper_1', 'f444w_bkg_aper_1', 'f444w_mask_aper_1', 'f444w_flux_aper_2', 'f444w_fluxerr_aper_2', 'f444w_flag_aper_2', 'f444w_bkg_aper_2', 'f444w_mask_aper_2', 'f444w_tot_corr', 'f606w_flux_aper_0', 'f606w_fluxerr_aper_0', 'f606w_flag_aper_0', 'f606w_bkg_aper_0', 'f606w_mask_aper_0', 'f606w_flux_aper_1', 'f606w_fluxerr_aper_1', 'f606w_flag_aper_1', 'f606w_bkg_aper_1', 'f606w_mask_aper_1', 'f606w_flux_aper_2', 'f606w_fluxerr_aper_2', 'f606w_flag_aper_2', 'f606w_bkg_aper_2', 'f606w_mask_aper_2', 'f606w_tot_corr', 'f606wu_flux_aper_0', 'f606wu_fluxerr_aper_0', 'f606wu_flag_aper_0', 'f606wu_bkg_aper_0', 'f606wu_mask_aper_0', 'f606wu_flux_aper_1', 'f606wu_fluxerr_aper_1', 'f606wu_flag_aper_1', 'f606wu_bkg_aper_1', 'f606wu_mask_aper_1', 'f606wu_flux_aper_2', 'f606wu_fluxerr_aper_2', 'f606wu_flag_aper_2', 'f606wu_bkg_aper_2', 'f606wu_mask_aper_2', 'f606wu_tot_corr', 'f814w_flux_aper_0', 'f814w_fluxerr_aper_0', 'f814w_flag_aper_0', 'f814w_bkg_aper_0', 'f814w_mask_aper_0', 'f814w_flux_aper_1', 'f814w_fluxerr_aper_1', 'f814w_flag_aper_1', 'f814w_bkg_aper_1', 'f814w_mask_aper_1', 'f814w_flux_aper_2', 'f814w_fluxerr_aper_2', 'f814w_flag_aper_2', 'f814w_bkg_aper_2', 'f814w_mask_aper_2', 'f814w_tot_corr', 'apcorr_0', 'f105w_corr_0', 'f105w_ecorr_0', 'f105w_tot_0', 'f105w_etot_0', 'f115w_corr_0', 'f115w_ecorr_0', 'f115w_tot_0', 'f115w_etot_0', 'f125w_corr_0', 'f125w_ecorr_0', 'f125w_tot_0', 'f125w_etot_0', 'f140w_corr_0', 'f140w_ecorr_0', 'f140w_tot_0', 'f140w_etot_0', 'f150w_corr_0', 'f150w_ecorr_0', 'f150w_tot_0', 'f150w_etot_0', 'f160w_corr_0', 'f160w_ecorr_0', 'f160w_tot_0', 'f160w_etot_0', 'f200w_corr_0', 'f200w_ecorr_0', 'f200w_tot_0', 'f200w_etot_0', 'f277w_corr_0', 'f277w_ecorr_0', 'f277w_tot_0', 'f277w_etot_0', 'f356w_corr_0', 'f356w_ecorr_0', 'f356w_tot_0', 'f356w_etot_0', 'f410m_corr_0', 'f410m_ecorr_0', 'f410m_tot_0', 'f410m_etot_0', 'f435w_corr_0', 'f435w_ecorr_0', 'f435w_tot_0', 'f435w_etot_0', 'f444w_corr_0', 'f444w_ecorr_0', 'f444w_tot_0', 'f444w_etot_0', 'f606w_corr_0', 'f606w_ecorr_0', 'f606w_tot_0', 'f606w_etot_0', 'f606wu_corr_0', 'f606wu_ecorr_0', 'f606wu_tot_0', 'f606wu_etot_0', 'f814w_corr_0', 'f814w_ecorr_0', 'f814w_tot_0', 'f814w_etot_0', 'apcorr_1', 'f105w_corr_1', 'f105w_ecorr_1', 'f105w_tot_1', 'f105w_etot_1', 'f115w_corr_1', 'f115w_ecorr_1', 'f115w_tot_1', 'f115w_etot_1', 'f125w_corr_1', 'f125w_ecorr_1', 'f125w_tot_1', 'f125w_etot_1', 'f140w_corr_1', 'f140w_ecorr_1', 'f140w_tot_1', 'f140w_etot_1', 'f150w_corr_1', 'f150w_ecorr_1', 'f150w_tot_1', 'f150w_etot_1', 'f160w_corr_1', 'f160w_ecorr_1', 'f160w_tot_1', 'f160w_etot_1', 'f200w_corr_1', 'f200w_ecorr_1', 'f200w_tot_1', 'f200w_etot_1', 'f277w_corr_1', 'f277w_ecorr_1', 'f277w_tot_1', 'f277w_etot_1', 'f356w_corr_1', 'f356w_ecorr_1', 'f356w_tot_1', 'f356w_etot_1', 'f410m_corr_1', 'f410m_ecorr_1', 'f410m_tot_1', 'f410m_etot_1', 'f435w_corr_1', 'f435w_ecorr_1', 'f435w_tot_1', 'f435w_etot_1', 'f444w_corr_1', 'f444w_ecorr_1', 'f444w_tot_1', 'f444w_etot_1', 'f606w_corr_1', 'f606w_ecorr_1', 'f606w_tot_1', 'f606w_etot_1', 'f606wu_corr_1', 'f606wu_ecorr_1', 'f606wu_tot_1', 'f606wu_etot_1', 'f814w_corr_1', 'f814w_ecorr_1', 'f814w_tot_1', 'f814w_etot_1', 'apcorr_2', 'f105w_corr_2', 'f105w_ecorr_2', 'f105w_tot_2', 'f105w_etot_2', 'f115w_corr_2', 'f115w_ecorr_2', 'f115w_tot_2', 'f115w_etot_2', 'f125w_corr_2', 'f125w_ecorr_2', 'f125w_tot_2', 'f125w_etot_2', 'f140w_corr_2', 'f140w_ecorr_2', 'f140w_tot_2', 'f140w_etot_2', 'f150w_corr_2', 'f150w_ecorr_2', 'f150w_tot_2', 'f150w_etot_2', 'f160w_corr_2', 'f160w_ecorr_2', 'f160w_tot_2', 'f160w_etot_2', 'f200w_corr_2', 'f200w_ecorr_2', 'f200w_tot_2', 'f200w_etot_2', 'f277w_corr_2', 'f277w_ecorr_2', 'f277w_tot_2', 'f277w_etot_2', 'f356w_corr_2', 'f356w_ecorr_2', 'f356w_tot_2', 'f356w_etot_2', 'f410m_corr_2', 'f410m_ecorr_2', 'f410m_tot_2', 'f410m_etot_2', 'f435w_corr_2', 'f435w_ecorr_2', 'f435w_tot_2', 'f435w_etot_2', 'f444w_corr_2', 'f444w_ecorr_2', 'f444w_tot_2', 'f444w_etot_2', 'f606w_corr_2', 'f606w_ecorr_2', 'f606w_tot_2', 'f606w_etot_2', 'f606wu_corr_2', 'f606wu_ecorr_2', 'f606wu_tot_2', 'f606wu_etot_2', 'f814w_corr_2', 'f814w_ecorr_2', 'f814w_tot_2', 'f814w_etot_2', 'z_spec', 'dummy_err', 'dummy_flux', 'nusefilt', 'z_ml', 'z_ml_chi2', 'z_ml_risk', 'lc_min', 'lc_max', 'z_phot', 'z_phot_chi2', 'z_phot_risk', 'z_min_risk', 'min_risk', 'z_raw_chi2', 'raw_chi2', 'z025', 'z160', 'z500', 'z840', 'z975', 'restU', 'restU_err', 'restB', 'restB_err', 'restV', 'restV_err', 'restJ', 'restJ_err', 'dL', 'mass', 'sfr', 'Lv', 'LIR', 'energy_abs', 'Lu', 'Lj', 'L1400', 'L2800', 'LHa', 'LOIII', 'LHb', 'LOII', 'MLv', 'Av', 'lwAgeV', 'rest120', 'rest120_err', 'rest121', 'rest121_err', 'rest156', 'rest156_err', 'rest157', 'rest157_err', 'rest158', 'rest158_err', 'rest159', 'rest159_err', 'rest160', 'rest160_err', 'rest414', 'rest414_err', 'rest415', 'rest415_err', 'rest416', 'rest416_err'
    ]


#_FLOAT_FEATURES = [
#    'thresh', 'npix', 'tnpix'
#    ]

class JWST(datasets.GeneratorBasedBuilder):
    """TODO: Short description of my dataset."""

    VERSION = _VERSION

    BUILDER_CONFIGS = [
        datasets.BuilderConfig(name="jwst_ceers_nircam", 
                               version=VERSION, 
                               data_files=DataFilesPatternsDict.from_patterns({'train': ['../../../../../data/astropile/ceers/ceers-full-grizli-v7.2_all_32_f115w-f150w-f200w-f277w-f356w-f444w/healpix=*/*.hdf5']}),
                               description="JWST NIRCam images CEERS sruvey"),
    ]

    DEFAULT_CONFIG_NAME = "jwst_ceers_nircam"

    _image_size = 96

    _bands = ['jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f200w','jwst_nircam_f277w','jwst_nircam_f356w','jwst_nircam_f444w']

    @classmethod
    def _info(self):
        """ Defines the features available in this dataset.
        """
        # Starting with all features common to image datasets
        features = {
            'image': Sequence(feature={
                'band': Value('string'),
                'array': Array2D(shape=(self._image_size, self._image_size), dtype='float32'),
                'psf_fwhm': Array2D(shape=(self._image_size, self._image_size), dtype='float32'),
                'scale': Value('float32'),
            })
        }
        # Adding all values from the catalog
        for f in _FLOAT_FEATURES:
            features[f] = Value('float32')

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
        data_files = dl_manager.download_and_extract(self.config.data_files)
        if isinstance(data_files, (str, list, tuple)):
            files = data_files
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            return [datasets.SplitGenerator(name=datasets.Split.TRAIN, gen_kwargs={"files": files})]
        splits = []
        for split_name, files in data_files.items():
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            splits.append(datasets.SplitGenerator(name=split_name, gen_kwargs={"files": files})) 
        return splits

    def _generate_examples(self, files, object_ids=None):
        """ Yields examples as (key, example) tuples.
        """
        for j, file in enumerate(itertools.chain.from_iterable(files)):
            print(file)
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                    
                else:
                    keys = data["object_id"]
                    #print(keys)
                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["object_id"])
                sorted_ids = data["object_id"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog 
                    i = sort_index[np.searchsorted(sorted_ids, k)]
                    # Parse image data
                    example = {'image':  [{'band': data['image_band'][i][j].decode('utf-8'),
                               'array': data['image_array'][i][j],
                               'psf_fwhm': data['image_psf_fwhm'][i][j],
                               'scale': data['image_scale'][i][j]} for j, _ in enumerate( self._bands )]
                    }
                    # Add all other requested features
                    for f in _FLOAT_FEATURES:
                        #print(f)
                        example[f] = data[f][i].astype('float32')
                    
                    # Add object_id
                    
                    example["object_id"] = str(data["object_id"][i])
                    #print(example)
                    
                    yield str(data['object_id'][i]), example
                    