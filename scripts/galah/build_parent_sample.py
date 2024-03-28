import os
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from multiprocessing import Pool
from tqdm import tqdm
import h5py
from collections import defaultdict
import pandas as pd
import healpy as hp

_PARAMETERS = [
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

# Following https://www.galah-survey.org/dr3/using_the_data/#recommended-flag-values
def selection_fn(catalog):
    mask = catalog['snr_c3_iraf'] > 30            # Reasonable signal-to-noise ratio
    mask &= catalog['flag_sp'] == 0
    mask &= catalog['flag_fe_h'] == 0
    mask &= catalog['survey_name'] == 'galah_main'
    return mask

def process_band_fits(filename):
    spectrum = {}
    try:
        hdus = fits.open(filename)
        
        # Unnormalized, sky-substracted spectrum
        flux = hdus[0].data
        sigma = hdus[1].data
        
        # Normalized spectrum
        norm_flux = hdus[4].data
        
        start_wavelength = hdus[0].header["CRVAL1"]
        dispersion       = hdus[0].header["CDELT1"]
        nr_pixels        = hdus[0].header["NAXIS1"]
        reference_pixel  = hdus[0].header["CRPIX1"]
        if reference_pixel == 0:
            reference_pixel = 1
        
        spectrum['flux'] = flux
        spectrum['lambda'] = ((np.arange(0,nr_pixels)--reference_pixel+1)*dispersion+start_wavelength)
        spectrum['ivar'] = 1 / np.power(flux * sigma, 2)
        
        # TODO: UPDATE
        # PLACEHOLDER VALUE!
        spectrum['lsf_sigma'] = np.ones_like(spectrum['lambda'])
        
        norm_start_wavelength = hdus[4].header["CRVAL1"]
        norm_dispersion       = hdus[4].header["CDELT1"]
        norm_nr_pixels        = hdus[4].header["NAXIS1"]
        norm_reference_pixel  = hdus[4].header["CRPIX1"]
        if norm_reference_pixel == 0:
            norm_reference_pixel=1
        spectrum['norm_flux'] = norm_flux
        spectrum['norm_lambda'] = ((np.arange(0,norm_nr_pixels)--norm_reference_pixel+1)*norm_dispersion+norm_start_wavelength)
        spectrum['norm_ivar'] = 1 / np.power(norm_flux * sigma, 2)
        
        return spectrum
    except FileNotFoundError as e:
        return defaultdict(list)

def processing_fn(args):
    filename_B, filename_G, filename_R, object_id = args
    
    spectrum_B, spectrum_G, spectrum_R = process_band_fits(filename_B), process_band_fits(filename_G), process_band_fits(filename_R)

    # Return the results
    return {'object_id': object_id,
            'spectrum_lambda': np.concatenate([spectrum_B['lambda'], spectrum_G['lambda'], spectrum_R['lambda']]),
            'spectrum_flux': np.concatenate([spectrum_B['flux'], spectrum_G['flux'], spectrum_R['flux']]),
            'spectrum_flux_ivar': np.concatenate([spectrum_B['ivar'], spectrum_G['ivar'], spectrum_R['ivar']]),
            'spectrum_lsf_sigma': np.concatenate([spectrum_B['lsf_sigma'], spectrum_G['lsf_sigma'], spectrum_R['lsf_sigma']]),
            'spectrum_norm_lambda': np.concatenate([spectrum_B['norm_lambda'], spectrum_G['norm_lambda'], spectrum_R['norm_lambda']]),
            'spectrum_norm_flux': np.concatenate([spectrum_B['norm_flux'], spectrum_G['norm_flux'], spectrum_R['norm_flux']]),
            'spectrum_norm_ivar': np.concatenate([spectrum_B['norm_ivar'], spectrum_G['norm_ivar'], spectrum_R['norm_ivar']])
            }


def save_in_standard_format(args):
    """ This function takes care of iterating through the different input files 
    corresponding to this healpix index, and exporting the data in standard format.
    """
    catalog, output_filename, galah_data_path = args
    catalog['object_id'] = catalog['sobject_id']
    
    # Create the output directory if it does not exist
    if not os.path.exists(os.path.dirname(output_filename)):
        os.makedirs(os.path.dirname(output_filename))

    # Process all files
    results = [
        processing_fn((
            os.path.join(galah_data_path, str(obj_id) + '1.fits'),
            os.path.join(galah_data_path, str(obj_id) + '2.fits'),
            os.path.join(galah_data_path, str(obj_id) + '3.fits'),
            obj_id
        )) for obj_id in catalog['object_id']
    ]
    
    # Pad all spectra to the same length
    max_length = max([len(d['spectrum_flux']) for d in results])
    for i in tqdm(range(len(results))):
        results[i]['spectrum_flux'] = np.pad(results[i]['spectrum_flux'], (0, max_length - len(results[i]['spectrum_flux'])), mode='constant')
        results[i]['spectrum_flux_ivar'] = np.pad(results[i]['spectrum_flux_ivar'], (0, max_length - len(results[i]['spectrum_flux_ivar'])), mode='constant')
        results[i]['spectrum_lambda'] = np.pad(results[i]['spectrum_lambda'], (0 , max_length - len(results[i]['spectrum_lambda'])), mode='constant', constant_values=-1)
        results[i]['spectrum_norm_flux'] = np.pad(results[i]['spectrum_norm_flux'], (0, max_length - len(results[i]['spectrum_norm_flux'])), mode='constant')
        results[i]['spectrum_norm_ivar'] = np.pad(results[i]['spectrum_norm_ivar'], (0, max_length - len(results[i]['spectrum_norm_ivar'])), mode='constant')
        results[i]['spectrum_norm_lambda'] = np.pad(results[i]['spectrum_norm_lambda'], (0, max_length - len(results[i]['spectrum_norm_lambda'])), mode='constant', constant_values=-1)
        results[i]['spectrum_flux'] = np.pad(results[i]['spectrum_flux'], (0, max_length - len(results[i]['spectrum_flux'])), mode='constant')

    # Aggregate all spectra into an astropy table
    spectra = Table({k: np.vstack([d[k] for d in results])
                     for k in results[0].keys()})
    spectra['object_id'] = spectra['object_id'].flatten()
    
    catalog = catalog[['object_id'] + _PARAMETERS]

    catalog = join(catalog, spectra, keys='object_id', join_type='right')

    # Making sure we didn't lose anyone
    assert len(catalog) == len(spectra), "There was an error in the join operation, probably some spectra files are missing"

    # Save all columns to disk in HDF5 format
    with h5py.File(output_filename, 'w') as hdf5_file:
        for key in catalog.colnames:
            hdf5_file.create_dataset(key, data=catalog[key])
    return 1

def main(args):
    # Load the catalog file and apply main cuts
    catalog = Table.read(os.path.join(args.galah_data_path, 'objects.csv'))
    
    catalog = catalog[selection_fn(catalog)]
    catalog['healpix'] = hp.ang2pix(64, catalog['ra'], catalog['dec'], lonlat=True, nest=True)
    
    h_catalog = catalog.group_by('healpix')
    
    map_args = []
    for group in h_catalog.groups:
        group_filename = os.path.join(args.output_dir, 'galah_dr3/healpix={}/001-of-001.hdf5'.format(group['healpix'][0]))
        map_args.append((group, group_filename, args.galah_data_path))
        
    with Pool(args.num_processes) as pool:
        results = list(tqdm(pool.imap(save_in_standard_format, map_args), total=len(map_args)))

    if sum(results) != len(map_args):
        print("There was an error in the parallel processing, some files may not have been processed correctly")

    print("All done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts spectra from all SDSS spectra downloaded through Globus')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    parser.add_argument('galah_data_path', type=str, help='Path to the local copy of the GALAH data')
    parser.add_argument('--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    args = parser.parse_args()

    main(args)