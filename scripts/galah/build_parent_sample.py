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
from scipy.optimize import curve_fit
from pathlib import Path
import requests
import tarfile
import requests
import pandas as pd
from typing import List
from waiting import wait

FILTERS = ["B", "G", "R", "I"]

def members(tar, strip):
    members = []
    object_ids = [f.path.split('/')[-1].split('.')[-2][:-1] for f in tar.getmembers() if f.path.endswith('.fits')][:2]
    for member in tar.getmembers():
        if member.path.endswith('.fits'):
            member.path = member.path.split('/', strip)[-1]
            if member.path.split('.')[0][:-1] in object_ids:
                members.append(member)
    return members

def get_object_ids(n_objects: int) -> pd.DataFrame:
    sql_query = """SELECT TOP {n_objects:d} * FROM galah_dr3.main_star WHERE sobject_id < 160000000000000 and logg < 2.0"""

    api_url = 'https://datacentral.org.au/api/services/query/'
    qdata = {
        'title' : 'galah_test_query',
        'sql' : sql_query.format(n_objects=n_objects),
        'run_async' : False
        }
    post = requests.post(api_url, data=qdata).json()
    resp = requests.get(post['url']).json()
    return pd.DataFrame(resp['result']['data'], columns=resp['result']['columns'])


def get_objects_table(source_ids: List[int] = None) -> pd.DataFrame:
    if source_ids:
        if not type(source_ids) == list and not all([type(si) == int for si in source_ids]):
            raise ValueError("source_ids must either be None or a list of integers")
        source_ids_str = "(" + ','.join([str(si) for si in source_ids]) + ")"
        sql_query = "SELECT * FROM galah_dr3.main_star WHERE sobject_id IN " + source_ids_str
    else:
        sql_query = "SELECT * FROM galah_dr3.main_star"
    api_url = 'https://datacentral.org.au/api/services/query/'
    qdata = {
        'title' : 'galah_test_query',
        'sql' : sql_query,
        'run_async' : True
        }
    post = requests.post(api_url,data=qdata).json()

    print("Downloading the GALAH DR3 objects table.")
    wait(lambda: requests.get(post['url']).json()['status']=='complete', timeout_seconds=600)
    resp = requests.get(post['url']).json()
    return pd.DataFrame(resp['result']['data'],columns=resp['result']['columns'])


def download_spectrum_file(object_id: int, filter_symbol: str) -> bytes:
    if filter_symbol not in FILTERS:
        raise Exception("Filter must be one of" + ",".join(FILTERS))
    
    url_template = "https://datacentral.org.au/vo/slink/links?ID={sobject_id:d}&DR=galah_dr3&IDX=0&FILT={filter}&RESPONSEFORMAT=fits"
    r = requests.get(url_template.format(sobject_id=object_id, filter=filter_symbol))
    return r.content


def download_resolution_maps() -> List[bytes]:
    url_template = "https://github.com/svenbuder/GALAH_DR3/raw/master/analysis/resolution_maps/ccd{ind:d}_piv.fits"
    return [requests.get(url_template.format(ind=i)).content for i in range(1, 5)]


FILTERS = ["B", "G", "R", "I"]

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

# These are going to be a 1d array in the example - one per example.
_FLOAT_FEATURES = [
    'object_id',
    'timestamp',
    'spectrum_B_ind_start',
    'spectrum_B_ind_end',
    'spectrum_G_ind_start',
    'spectrum_G_ind_end',
    'spectrum_R_ind_start',
    'spectrum_R_ind_end',
    'spectrum_I_ind_start',
    'spectrum_I_ind_end'
]

# Following https://www.galah-survey.org/dr3/using_the_data/#recommended-flag-values
def selection_fn(catalog):
    mask = catalog['snr_c3_iraf'] > 30            # Reasonable signal-to-noise ratio
    mask &= catalog['flag_sp'] == 0
    mask &= catalog['flag_fe_h'] == 0
    mask &= catalog['survey_name'] == 'galah_main'
    return mask


def get_resolution(ccd_resolution_map_filename):
    try:
        resolution = fits.open(ccd_resolution_map_filename)[0]
    except IndexError:
        raise IOError("The resolution map file is invalid or empty!")
    
    mean_resolution = np.mean(resolution.data, axis=0)
    def _gauss(x, a, x0, sigma):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    try:
        popt, _ = curve_fit(_gauss, np.arange(len(mean_resolution)), mean_resolution, p0=[1, 5, 1])
    except:
        popt = [-1, -1, -1]
    return mean_resolution, popt[2]*np.ones(shape=len(mean_resolution), dtype=np.float32)


def process_band_fits(filename, resolution_filename):
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
            # Get an averaged estimated Gaussian line spread function
        # TODO: Actually properly estimate the line spread function of each spectrum
        lsf, lsf_sigma = get_resolution(resolution_filename)
        
        norm_start_wavelength = hdus[4].header["CRVAL1"]
        norm_dispersion       = hdus[4].header["CDELT1"]
        norm_nr_pixels        = hdus[4].header["NAXIS1"]
        norm_reference_pixel  = hdus[4].header["CRPIX1"]
        if norm_reference_pixel == 0:
            norm_reference_pixel=1
        spectrum['norm_flux'] = norm_flux
        spectrum['norm_lambda'] = ((np.arange(0,norm_nr_pixels)--norm_reference_pixel+1)*norm_dispersion+norm_start_wavelength)
        spectrum['norm_ivar'] = 1 / np.power(norm_flux * sigma, 2)
        spectrum['lsf'] = lsf
        spectrum['lsf_sigma'] = lsf_sigma
        spectrum['timestamp'] = hdus[0].header['UTMJD']
        
        return spectrum
    except FileNotFoundError as e:
        return defaultdict(list)
    
def download_galah_data(output_path, tiny=False):
    if not Path(output_path).exists():
        Path(output_path).mkdir(parents=True)
    
    if tiny:
        # Download a small sample of GALAH data for testing
        catalog_url = "https://cloud.datacentral.org.au/teamdata/GALAH/public/GALAH_DR3/spectra/tar_files/131116_com.tar.gz"
        response = requests.get(catalog_url, stream=True)

        # Sizes in bytes.
        total_size = int(response.headers.get("content-length", 0))
        block_size = 1024

        with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
            with open(os.path.join(output_path, 'spectra.tar.gz'), "wb") as file:
                for data in response.iter_content(block_size):
                    progress_bar.update(len(data))
                    file.write(data)

        if total_size != 0 and progress_bar.n != total_size:
            raise RuntimeError("Could not download file")
        
        zipped_file = tarfile.open(os.path.join(output_path, 'spectra.tar.gz'))
        strip = max([len(n.split('/')) for n in zipped_file.getnames() if n.endswith('.fits')])
        zipped_file.extractall(output_path, members=members(zipped_file, strip=strip)[:8])
        zipped_file.close()
        os.remove(zipped_file.name)
        
        filenames = os.listdir(output_path)
        # move the spectra
        object_ids = list(set([int(f.split('.')[0][:-1]) for f in filenames if '.fits' in f]))

        objects_data = get_objects_table(object_ids)
        objects_data.to_csv(os.path.join(output_path, 'objects.csv'), index=False)
    
    else:
        object_data = get_objects_table()
        object_data.to_csv(os.path.join(output_path, 'objects.csv'), index=False)
        
        catalog_url = "https://cloud.datacentral.org.au/teamdata/GALAH/public/GALAH_DR3/spectra/GALAH_DR3_all_spectra_with_normalisation_v2.tar.gz"
        # Streaming, so we can iterate over the response.
        response = requests.get(catalog_url, stream=True)

        # Sizes in bytes.
        total_size = int(response.headers.get("content-length", 0))
        block_size = 1024

        with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
            with open(os.path.join(output_path, 'spectra.tar.gz'), "wb") as file:
                for data in response.iter_content(block_size):
                    progress_bar.update(len(data))
                    file.write(data)

        if total_size != 0 and progress_bar.n != total_size:
            raise RuntimeError("Could not download file")
        
        zipped_file = tarfile.open(os.path.join(output_path, 'spectra.tar.gz'))
        zipped_file.extractall(output_path)
        zipped_file.close()
        os.remove(zipped_file.name)
        
    # Resolution maps
    resolution_maps_path = Path(os.path.join(output_path, "resolution_maps"))
    resolution_maps_path.mkdir(parents=True)
    
    resolution_maps = download_resolution_maps()
    for i, rm in tqdm(enumerate(resolution_maps)):
        open(os.path.join(resolution_maps_path, "ccd{i:d}_piv.fits".format(i=i+1)), "wb").write(rm)
    

def processing_fn(args):
    (filename_B, filename_G, filename_R, filename_I,
     resolution_B, resolution_G, resolution_R, resolution_I,
     object_id) = args
    
    spectrum_B, spectrum_G, spectrum_R, spectrum_I = process_band_fits(filename_B, resolution_B), \
        process_band_fits(filename_G, resolution_G), process_band_fits(filename_R, resolution_R), \
        process_band_fits(filename_I, resolution_I)
        
    len_B = len(spectrum_B['lambda'])
    len_G = len(spectrum_G['lambda'])
    len_R = len(spectrum_R['lambda'])
    len_I = len(spectrum_I['lambda'])
    
    def concatenate_spectra(spectra_key: str):
        return np.concatenate([spectrum_B[spectra_key], spectrum_G[spectra_key], spectrum_R[spectra_key], spectrum_I[spectra_key]]).flatten()

    # Return the results
    return {'object_id': object_id,
            'timestamp': np.mean([spectrum_B['timestamp'], spectrum_G['timestamp'], spectrum_R['timestamp'], spectrum_I['timestamp']]),
            'spectrum_lambda': concatenate_spectra('lambda'),
            'spectrum_flux': concatenate_spectra('flux'),
            'spectrum_flux_ivar': concatenate_spectra('ivar'),
            'spectrum_lsf': concatenate_spectra('lsf'),
            'spectrum_lsf_sigma': concatenate_spectra('lsf_sigma'),
            'spectrum_norm_lambda': concatenate_spectra('norm_lambda'),
            'spectrum_norm_flux': concatenate_spectra('norm_flux'),
            'spectrum_norm_ivar': concatenate_spectra('norm_ivar'),
            'spectrum_B_ind_start': 0,
            'spectrum_B_ind_end': len_B-1,
            'spectrum_G_ind_start': len_B,
            'spectrum_G_ind_end': len_B+len_G-1,
            'spectrum_R_ind_start': len_B+len_G,
            'spectrum_R_ind_end': len_B+len_G+len_R-1,
            'spectrum_I_ind_start': len_B+len_G+len_R,
            'spectrum_I_ind_end': len_B+len_G+len_R+len_I-1
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
            os.path.join(galah_data_path, str(obj_id) + '4.fits'),
            os.path.join(galah_data_path, 'resolution_maps/ccd1_piv.fits'),
            os.path.join(galah_data_path, 'resolution_maps/ccd2_piv.fits'),
            os.path.join(galah_data_path, 'resolution_maps/ccd3_piv.fits'),
            os.path.join(galah_data_path, 'resolution_maps/ccd4_piv.fits'),
            obj_id
        )) for obj_id in catalog['object_id']
    ]
    
    # Pad all spectra to the same length
    max_length = max([len(d['spectrum_flux']) for d in results])
    for i in tqdm(range(len(results))):
        results[i]['spectrum_flux'] = np.pad(results[i]['spectrum_flux'], (0, max_length - len(results[i]['spectrum_flux'])), mode='constant')
        results[i]['spectrum_flux_ivar'] = np.pad(results[i]['spectrum_flux_ivar'], (0, max_length - len(results[i]['spectrum_flux_ivar'])), mode='constant')
        results[i]['spectrum_lsf'] = np.pad(results[i]['spectrum_lsf'], (0, max_length - len(results[i]['spectrum_lsf'])), mode='constant')
        results[i]['spectrum_lsf_sigma'] = np.pad(results[i]['spectrum_lsf_sigma'], (0, max_length - len(results[i]['spectrum_lsf_sigma'])), mode='constant')
        results[i]['spectrum_lambda'] = np.pad(results[i]['spectrum_lambda'], (0 , max_length - len(results[i]['spectrum_lambda'])), mode='constant', constant_values=-1)
        results[i]['spectrum_norm_flux'] = np.pad(results[i]['spectrum_norm_flux'], (0, max_length - len(results[i]['spectrum_norm_flux'])), mode='constant')
        results[i]['spectrum_norm_ivar'] = np.pad(results[i]['spectrum_norm_ivar'], (0, max_length - len(results[i]['spectrum_norm_ivar'])), mode='constant')
        results[i]['spectrum_norm_lambda'] = np.pad(results[i]['spectrum_norm_lambda'], (0, max_length - len(results[i]['spectrum_norm_lambda'])), mode='constant', constant_values=-1)

    # Aggregate all spectra into an astropy table
    spectra = Table({k: np.vstack([d[k] for d in results])
                     for k in results[0].keys()})
    
    for feature in _FLOAT_FEATURES:
        spectra[feature] = spectra[feature].flatten()
    
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
    if not os.path.exists(args.galah_data_path):
        Path(args.galah_data_path).mkdir(parents=True)

    if not os.listdir(args.galah_data_path):
        download_galah_data(args.galah_data_path, args.tiny)
    
    catalog = Table.read(os.path.join(args.galah_data_path, 'objects.csv'))
    
    if not args.tiny:
        catalog = catalog[selection_fn(catalog)]
    catalog['healpix'] = hp.ang2pix(64, catalog['ra'], catalog['dec'], lonlat=True, nest=True)
    
    h_catalog = catalog.group_by('healpix')
    
    map_args = []
    for group in h_catalog.groups:
        group_filename = os.path.join(args.output_path, 'galah_dr3/healpix={}/001-of-001.hdf5'.format(group['healpix'][0]))
        map_args.append((group, group_filename, args.galah_data_path))
        
    with Pool(args.num_processes) as pool:
        results = list(tqdm(pool.imap(save_in_standard_format, map_args), total=len(map_args)))

    if sum(results) != len(map_args):
        print("There was an error in the parallel processing, some files may not have been processed correctly")

    print("All done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts spectra from GALAH')
    parser.add_argument('galah_data_path', type=str, help='Path to the local copy of the GALAH data')
    parser.add_argument('output_path', type=str, help='Path to the output directory')
    parser.add_argument('--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    parser.add_argument('--tiny', action='store_true', help='Use a tiny subset of the data for testing')
    args = parser.parse_args()

    main(args)