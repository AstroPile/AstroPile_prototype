
import argparse
import itertools
import pathlib
from multiprocessing import Pool
import numpy as np

from astropy.io import fits
from astropy.table import Table, join
from tqdm import tqdm
import healpy as hp
import h5py


_utf8_filter_type = h5py.string_dtype('utf-8', 5)
_healpix_nside = 16

def get_maps_data(plateifu: str, cubefile: pathlib.Path, nspaxels: int, pad_arr: tuple) -> dict:
    """ Extract MaNGA DAP MAPs data

    Extract information from the MaNGA DAP MAPS-HYB10-MILESHC-MASTARSSP files.
    Gets the sky x,y coords; elliptical r,theta coordinates, and r/re and r_kpc radii.

    Parameters
    ----------
    plateifu : str
        the MaNGA plateifu
    cubefile : pathlib.Path
        the DRP LOGCUBE filepath
    nspaxels : int
        the number of spaxels, this already padded to 96
    pad_arr : tuple
        the array of padding values

    Returns
    -------
    dict
        the map data
    """
    daptype = 'HYB10-MILESHC-MASTARSSP'
    base = cubefile.parents[5] / 'spectro/analysis/v3_1_1/3.1.0'
    plate, ifu = plateifu.split('-')
    mapfile = base / daptype / plate / ifu / f'manga-{plateifu}-MAPS-{daptype}.fits.gz'

    # open the DAP MAPS-HYB10-MILESHC-MASTARSSP file
    with fits.open(mapfile) as hdulist:
        data = {'coords': {}}

        # get spaxel coordinates; pad image shape to 96
        skycoo_x, skycoo_y = np.pad(hdulist['SPX_SKYCOO'].data, pad_arr).astype(np.float32)
        skycoo_units = hdulist['SPX_SKYCOO'].header['BUNIT'].encode('utf-8')
        ellcoo_r, ellcoo_rre, ellcoo_rkpc, ellcoo_theta = np.pad(hdulist['SPX_ELLCOO'].data, pad_arr).astype(np.float32)
        ellcoo_r_units, ellcoo_rre_units, ellcoo_rkpc_units, ellcoo_theta_units = \
            [hdulist['SPX_ELLCOO'].header.get(f'U{i}').encode('utf-8') for i in range(1, 5)]

        # reshape coord arrs to nspaxels
        data["coords"]['skycoo_x'] = skycoo_x.reshape(nspaxels)
        data["coords"]['skycoo_y'] = skycoo_y.reshape(nspaxels)
        data["coords"]['ellcoo_r'] = ellcoo_r.reshape(nspaxels)
        data["coords"]['ellcoo_rre'] = ellcoo_rre.reshape(nspaxels)
        data["coords"]['ellcoo_rkpc'] = ellcoo_rkpc.reshape(nspaxels)
        data["coords"]['ellcoo_theta'] = ellcoo_theta.reshape(nspaxels)
        data["coords"]['skycoo_units'] = itertools.repeat(skycoo_units)
        data["coords"]['ellcoo_r_units'] = itertools.repeat(ellcoo_r_units)
        data["coords"]['ellcoo_rre_units'] = itertools.repeat(ellcoo_rre_units)
        data["coords"]['ellcoo_rkpc_units'] = itertools.repeat(ellcoo_rkpc_units)
        data["coords"]['ellcoo_theta_units'] = itertools.repeat(ellcoo_theta_units)

        # grab map data
        maps = []
        for ext in hdulist:
            if ext.name == 'PRIMARY' or ext.name.endswith(("IVAR", "MASK")):
                continue

            name = ext.name.lower()
            unit = ext.header.get("BUNIT", '')
            array = ext.data
            ndim = ext.data.ndim

            # set padding array and channels based on ndim
            if ndim == 3:
                parr = pad_arr
                nchannels = array.shape[0]
            elif ndim == 2:
                parr = pad_arr[1:]
                nchannels = None

            # pad the data
            padded_data = np.pad(array, parr)
            dshape = padded_data.shape

            # get the uncertainty and mask arrays if any
            # pad images to 96 shape
            # pad mask with 1073741824 values as DONOTUSE
            errdata = ext.header.get("ERRDATA")
            qualdata = ext.header.get("QUALDATA")
            padded_errdata = np.pad(hdulist[errdata].data, parr) if errdata else np.zeros(dshape)
            padded_qualdata = np.pad(hdulist[qualdata].data, parr, constant_values=1073741824) if qualdata else np.full(dshape, 1073741824)

            # grab the maps, expand individual channels
            if nchannels:
                for channel in range(nchannels):
                    # format channel and units names
                    channame = ext.header.get(f'C{channel + 1:02}', ext.header.get(f'C{channel + 1}', '')).replace('-', '_').strip().replace('. ', '').replace(' ', '_')
                    uname = ext.header.get(f'U{channel+1:02}', ext.header.get(f'U{channel+1}', '')).replace('-', '_').strip().replace('. ', '').replace(' ', '_')
                    label = f'{ext.name.lower()}_{channame.lower()}'
                    unit = uname or unit

                    maps.append({
                        "group": name.encode('utf-8'),
                        "label": label.encode('utf-8'),
                        "array": padded_data[channel, :].astype(np.float32),
                        "ivar": padded_errdata[channel, :].astype(np.float32),
                        "mask": padded_qualdata[channel, :].astype(np.float32),
                        "array_units": unit.encode('utf-8')
                    })

            else:
                maps.append({
                    "group": name.encode('utf-8'),
                    "label": name.encode('utf-8'),
                    "array": padded_data.astype(np.float32),
                    "ivar": padded_errdata.astype(np.float32),
                    "mask": padded_qualdata.astype(np.float32),
                    "array_units": unit.encode('utf-8')
                })

        data['maps'] = maps

        return data


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

        # create x, y array indices and unique spaxel index
        y, x = np.indices((nx, ny))
        x = x.reshape(1, nspaxels)
        y = y.reshape(1, nspaxels)
        spaxel_idx = np.arange(nspaxels)

        # conversion from spaxel_idx back to x, y
        # y = int(spaxel_idx / float(nx))
        # x = spaxel_idx - y * nx

        # reshape and grab arrays
        # pad mask array with 1024 to indicate as DONOTUSE
        # pad rest with 0s
        flux = flux.reshape(nwave, nspaxels)
        ivar = np.pad(hdulist['IVAR'].data, pad_arr).reshape(nwave, nspaxels)
        mask = np.pad(hdulist['MASK'].data, pad_arr, constant_values=1024).reshape(nwave, nspaxels)
        lsf = np.pad(hdulist['LSFPOST'].data, pad_arr).reshape(nwave, nspaxels)

        wave = hdulist['WAVE'].data.astype(np.float32)
        wave = np.repeat(wave[:, np.newaxis], nspaxels, axis=1)

        # get DAP map data
        mapdata = get_maps_data(summary_row['plateifu'], filename, nspaxels, pad_arr)

        # add spaxels
        # combine the spaxel data together
        keys = ['flux', 'ivar', 'mask', 'lsf_sigma', 'lambda', 'x', 'y', 'spaxel_idx', 'flux_units', 'lambda_units']
        zz = zip(flux.T, ivar.T, mask.T, lsf.T, wave.T, x[0, :], y[0, :], spaxel_idx, flux_units, lambda_units)

        # optionally add any mapdata
        if mapdata:
            keys.extend(mapdata['coords'].keys())
            bb = zip(*mapdata['coords'].values())
            zz = zip(*(list(zip(*zz)) + list(zip(*bb))))

        # convert spaxels to a list of dicts
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

        # add the DAP maps data
        data['maps'] = mapdata['maps']

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
        plate, _ = plateifu.split('-')
        file = pathlib.Path(data_path) / 'dr17/manga/spectro/redux/v3_1_1' / plate / 'stack' / f'manga-{plateifu}-LOGCUBE.fits.gz'
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
                 if key not in ('provenance', 'spaxels', 'images', 'maps'):
                     hg.attrs[key] = res[key]
                     hg.create_dataset(key, data=res[key])

            # load spaxels
            spax = Table(res['spaxels'])
            hg.create_dataset('spaxels', data=spax)

            # load images
            im = Table(res['images'][0])
            hg.create_dataset('images', data=im)

            # load the maps data
            maps = Table(res['maps'])
            hg.create_dataset('maps', data=maps)

    return 1


def process_files(manga_data_path: str, output_dir: str, num_processes: int = 10, tiny=False):
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
    catalog = Table.read(manga_data_path + '/' + 'drpall-v3_1_1.fits', hdu='MANGA')
    catalog_dap = Table.read(manga_data_path + '/' + 'dapall-v3_1_1-3.1.0.fits', hdu='HYB10-MILESHC-MASTARSSP')
    catalog = join(catalog, catalog_dap, keys_left='plateifu', keys_right='PLATEIFU',  join_type='inner')
    # Removing entries for which the DAP was not completed
    catalog = catalog[catalog['DAPDONE']]

    if tiny:
        m = catalog['plateifu'] == '8485-1901'
        catalog = catalog[m]

    # Add healpix index to the catalog, and group the table
    catalog['healpix'] = hp.ang2pix(_healpix_nside, catalog['ifura'], catalog['ifudec'], lonlat=True, nest=True)
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

    # if sum(results) != len(map_args):
    #     print("There was an error in the parallel processing, some files may not have been processed correctly")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts IFU data from all SDSS MaNGA downloaded')
    parser.add_argument('-m', '--manga_data_path', type=str, help='Path to the local copy of the SDSS MaNGA data')
    parser.add_argument('-o', '--output_dir', type=str, default='out', help='Path to the output directory')
    parser.add_argument('-n', '--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    parser.add_argument('--tiny', action="store_true", help='Use a small subset of the data for testing')
    args = parser.parse_args()

    process_files(args.manga_data_path, args.output_dir, args.num_processes, args.tiny)