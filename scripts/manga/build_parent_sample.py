
import argparse
import pathlib
from multiprocessing import Pool
import numpy as np

from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
import healpy as hp
import h5py


_utf8_filter_type = h5py.string_dtype('utf-8', 5)

# def selection_fn(catalog):
#     return catalog


def process_single_plateifu(args: tuple) -> dict:
    """ Process a single MaNGA plate-IFU

    Extract the relevant information for a SDSS MaNGA
    plate-IFU observation.  Pulls from the MaNGA LOGCUBE fits
    file.  Extracts the plateifu as observation id, RA, Dec, healpix
    id, redshift, and spaxel size.  Also extracts the spaxel and
    reconstructed image data.  We pad the IFU data up to the image size
    of 96 elements.  Typical cube sizes in MaNGA range from 20-80 array
    elements.  This paddinng also adds empty spaxels.

    For each spaxel, we include the flux, ivar, lsf, wavelength, x and y
    array indices, and flux/wave units.

    For image data, we include the reconstructed filter-band and PSF images,
    the filter band, and image pixel units.

    Parameters
    ----------
    args : tuple
        input arguments

    Returns
    -------
    dict
        the output extracted data
    """

    summary_row, filename, object_id = args

    # set up data object
    data = {}
    data["provenance"] = {"project": "SDSS", "survey": "MaNGA", "release": "DR17"}

    # add meta data
    data['object_id'] = summary_row['plateifu']
    data['ra'] = summary_row['ifura']
    data['dec'] = summary_row['ifudec']
    data['healpix'] = summary_row['healpix']
    data['z'] = summary_row['nsa_z']
    data['spaxel_size'] = 0.5
    data['spaxel_size_unit'] = b'arcsec'

    # Load the CUBE file
    with fits.open(filename) as hdulist:
        flux = hdulist['FLUX'].data
        nwave, ny, nx = flux.shape
        nspaxels = ny * nx

        # compute padding size ; pad up to 96
        # pad_arr is 0 for z-axis, and padding for (before, after) in each spatial axis
        padding = int((96 - nx) / 2)
        pad_arr = ((0, 0), (padding, padding), (padding, padding))

        # repad flux
        flux = np.pad(flux, pad_arr)
        nwave, ny, nx = flux.shape
        nspaxels = ny * nx

        # units
        flux_units = hdulist['FLUX'].header['BUNIT'].encode('utf-8')
        lambda_units = hdulist['FLUX'].header['CUNIT3'].encode('utf-8')
        flux_units = np.repeat(flux_units, nspaxels)
        lambda_units = np.repeat(lambda_units, nspaxels)

        # create x, y array indices
        y, x = np.indices((nx, ny))
        x = x.reshape(1, nspaxels)
        y = y.reshape(1, nspaxels)

        # reshape and grab arrays
        # pad mask array with 1024 to indicate as DONOTUSE
        # pad rest with 0s
        flux = flux.reshape(nwave, nspaxels)
        ivar = np.pad(hdulist['IVAR'].data, pad_arr).reshape(nwave, nspaxels)
        mask = np.pad(hdulist['FLUX'].data, pad_arr, constant_values=1024).reshape(nwave, nspaxels)
        lsf = np.pad(hdulist['LSFPOST'].data, pad_arr).reshape(nwave, nspaxels)

        wave = hdulist['WAVE'].data.astype(np.float32)
        wave = np.repeat(wave[:, np.newaxis], nspaxels, axis=1)

        # add spaxels

        # combine the data together
        keys = ['flux', 'ivar', 'mask', 'lsf_sigma', 'lambda', 'x', 'y', 'flux_units', 'lambda_units']
        zz = zip(flux.T, ivar.T, mask.T, lsf.T, wave.T, x[0, :], y[0, :], flux_units, lambda_units)
        spaxels = [dict(zip(keys, values)) for values in zz]
        data['spaxels'] = spaxels

        # add images
        images = []
        n_filters = 4
        filters = np.array(['g', 'r', 'i', 'z'], dtype=_utf8_filter_type)
        img = create_images(hdulist, 'img', pad_arr[1:])
        psf = create_images(hdulist, 'psf', pad_arr[1:])
        images.append({
            "image_band": filters,
            "image_array": img,
            "image_array_units": np.array([b"nanomaggies/pixel"] * n_filters),
            "image_psf": psf,
            "image_psf_units": np.array([b"nanomaggies/pixel"] * n_filters),
            "image_scale": np.array([0.5] * n_filters),
            "image_scale_units": np.array([b"arcsec"] * n_filters)
        })
        data['images'] = images

    # Return the results
    return data


def create_images(hdu: fits.HDUList, image_type: str, pad_arr: tuple) -> np.array:
    """ Create a stack of images from the MaNGA data

    From the given MaNGA plate-IFU observation, extracts the reconstructed
    image data and stacks all the filters (g, r, i, z) together into a single
    array. Image type can be "img" to extracted the reconstructed filter image
    or "psf" to extract the reconstructed PSF in that filter.

    Parameters
    ----------
    hdu : fits.HDUList
        the cube fits data for the plate-IFU
    image_type : str
        the type of image data to extract
    pad_arr : tuple
        the padding size for each dimension

    Returns
    -------
    np.array
        the stacked image array
    """
    if image_type == 'img':
        g = np.pad(hdu['GIMG'].data, pad_arr).astype(np.float32)
        r = np.pad(hdu['RIMG'].data, pad_arr).astype(np.float32)
        i = np.pad(hdu['IIMG'].data, pad_arr).astype(np.float32)
        z = np.pad(hdu['ZIMG'].data, pad_arr).astype(np.float32)
    elif image_type == 'psf':
        g = np.pad(hdu['GPSF'].data, pad_arr).astype(np.float32)
        r = np.pad(hdu['RPSF'].data, pad_arr).astype(np.float32)
        i = np.pad(hdu['IPSF'].data, pad_arr).astype(np.float32)
        z = np.pad(hdu['ZPSF'].data, pad_arr).astype(np.float32)

    img_arr = np.stack([g, r, i, z])
    return img_arr


def process_healpix_group(args: tuple) -> int:
    """ Process a healpix group

    Process a group of plate-IFUS by healpix id.  The input args
    is a tuple of the (astropy.Table group, the output hdf5 filename
    for the group, and the input data path).  Writes the processed
    plate-IFUS into the designated HDF5 file.

    Parameters
    ----------
    args : tuple
        input arguments

    Returns
    -------
    int
        1 or 0 for success or failure
    """
    hp_grp, output_filename, data_path = args

    # Create the output directory if it does not exist
    path = pathlib.Path(output_filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Preparing the arguments for the parallel processing
    map_args = []
    for row in hp_grp:
        plateifu = row['plateifu']
        # find the data cube filepath
        files = pathlib.Path(data_path).rglob(f'*{plateifu}*LOGCUBE.fits*')
        for file in files:
            if file.exists():
                map_args.append((row, file, plateifu))

    # Process all files
    results = []
    for args in map_args:
        results.append(process_single_plateifu(args))

    if not results:
        return 0

    # Save all results to disk in HDF5 format
    with h5py.File(output_filename, 'w') as hdf:
        prov = results[0]['provenance']
        hdf.attrs['project'] = prov['project']
        hdf.attrs['survey'] = prov['survey']
        hdf.attrs['release'] = prov['release']

        for res in results:
            obsid = res['object_id']
            hdf.create_group(obsid, track_order=True)
            hg = hdf[obsid]

            # load metadata
            for key in res.keys():
                 if key not in ('provenance', 'spaxels', 'images'):
                     hg.attrs[key] = res[key]
                     hg.create_dataset(key, data=res[key])

            # load spaxels
            spax = Table(res['spaxels'])
            hg.create_dataset('spaxels', data=spax)

            # load images
            #im = Table({k: [d[k] for d in res['images']] for k in res['images'][0].keys()})
            im = Table(res['images'][0])
            hg.create_dataset('images', data=im)

    return 1


def process_files(manga_data_path: str, output_dir: str, num_processes: int = 10):
    """ Process SDSS MaNGA files

    Process downloaded SDSS MaNGA files using multiprocessing parallelization.
    Organizes the MaNGA drpall catalog by healpix id and processes plate-IFUs
    by healpix groups.  Within the output_dir path, files are organized in
    directories manga/healpix=****/001-of-001.hdf5.

    Parameters
    ----------
    manga_data_path : str
        the top level directory to the data
    output_dir : str
        the output directory for the hdf5 files
    num_processes : int, optional
        the number of processess to use, by default 10
    """
    # Load the catalog file and apply main cuts
    path = pathlib.Path(manga_data_path) / 'drpall-v3_1_1.fits'
    catalog = Table.read(path, hdu='MANGA')
    #catalog = catalog[selection_fn(catalog)]

    # Add healpix index to the catalog, and group the table
    catalog['healpix'] = hp.ang2pix(64, catalog['ifura'], catalog['ifudec'], lonlat=True, nest=True)
    hp_groups = catalog.group_by(['healpix'])

    # Preparing the arguments for the parallel processing
    map_args = []
    for group in hp_groups.groups:
        # Create a filename for the group
        path = pathlib.Path(output_dir) / f'manga/healpix={group["healpix"][0]}/001-of-001.hdf5'
        map_args.append((group, path, manga_data_path))

    # Run the parallel processing
    with Pool(num_processes) as pool:
        results = list(tqdm(pool.imap(process_healpix_group, map_args), total=len(map_args)))

    if sum(results) != len(map_args):
        print("There was an error in the parallel processing, some files may not have been processed correctly")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts IFU data from all SDSS MaNGA downloaded')
    parser.add_argument('manga_data_path', type=str, help='Path to the local copy of the SDSS MaNGA data')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    parser.add_argument('--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    args = parser.parse_args()

    process_files(args)