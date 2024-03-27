import os
import pathlib
import argparse
import glob
import numpy as np
from astropy.io import fits
from multiprocessing import Pool
import tarfile
from tqdm import tqdm
from urllib.request import urlretrieve
import h5py

def download_data(vipers_data_path):
    print ("Downloading files (unless available locally)")

    # Create the output directory if it does not exist
    if not os.path.exists(vipers_data_path):
        os.makedirs(vipers_data_path)

    url = "http://vipers.inaf.it/data/pdr2/spectra/"
    files = ["VIPERS_W1_SPECTRA_1D_PDR2.tar.gz", "VIPERS_W4_SPECTRA_1D_PDR2.tar.gz"]

    for f in files:
        local_path = os.path.join(vipers_data_path, f)
        if not os.path.exists(local_path):
            urlretrieve(url + f, local_path)

        tar = tarfile.open(local_path)
        tar.extractall(filter='data', path=vipers_data_path)
        tar.close()

def extract_data(filename):
    hdu = fits.open(filename)
    header = hdu[1].header
    id, ra, dec, redshift = int(header['ID']), header['RA'], header['DEC'], header['REDSHIFT']
    # TODO: has also redshift flag, exptime, norm, and mag. Needed?
    spec = hdu[1].data['FLUXES']
    wave = hdu[1].data['WAVES']
    noise = hdu[1].data['NOISE']
    mask = hdu[1].data['MASK']
    hdu.close()
    return id, ra, dec, redshift, spec, wave, noise, mask

def save_to_file(results, filename):
    print("Saving HDF5 file")
    keys = ["ID", "RA", "DEC", "REDSHIFT", "WAVELENGTH", "FLUX", "NOISE", "MASK"]
    results = {
        key: np.stack([d[i] for d in results], axis=0)
        for i, key in enumerate(keys)
    }

    with h5py.File(filename, "w") as hdf5_file:
        for key in keys:
            hdf5_file.create_dataset(key, data=results[key])


def main(args):

    # Download the data
    download_data(args.vipers_data_path)

    # Load all fits file, standardize them and append to HDF5 file
    files = glob.glob(os.path.join(args.vipers_data_path, '*.fits'))
    files = files

    # Run the parallel processing
    with Pool(args.num_processes) as pool:
        results = list(tqdm(pool.imap(extract_data, files), total=len(files)))

    if len(results) != len(files):
        print("There was an error in the parallel processing, some files may not have been processed correctly")

    # Create HDF5 file in standard format
    save_to_file(results, os.path.join(args.output_dir, 'vipers_dr2.hdf5'))

    # Remove all temp fits files, keep HDF5 and tar.gz
    print ("Cleaning up")
    for f in files:
        pathlib.Path(f).unlink(missing_ok=True)

    print("All done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts spectra from all VIPERS spectra downloaded from the web')
    parser.add_argument('vipers_data_path', type=str, help='Path to the local copy of the VIPERS data')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    parser.add_argument('--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    args = parser.parse_args()

    main(args)