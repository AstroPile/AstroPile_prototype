import os
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, join, hstack
from multiprocessing import Pool
from tqdm import tqdm
import healpy as hp
import h5py
import urllib

# APOGEE shares a global wavelength grid
lam = 10.0 ** np.arange(
    4.179, 4.179 + 8575 * 6.0 * 10.0**-6.0, 6.0 * 10.0**-6.0
).astype(np.float32)

# detector edges
blue_start = 246
blue_end = 3274
green_start = 3585
green_end = 6080
red_start = 6344
red_end = 8335
lam_cropped = lam[np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]]


def selection_fn(base_path, catalog):
    # Only use the spectrum from APO 2.5m and LCO 2.5m
    mask = (catalog["TELESCOPE"] == "apo25m") | (catalog["TELESCOPE"] == "lco25m")
    # known no file entries
    mask &= ~catalog["FILE"].mask

    # go thru all to check if those files actually exist (some files are missing even on APOGEE server)
    for idx, (apogee_id, field, telescope, filename) in enumerate(zip(catalog["APOGEE_ID"], catalog["FIELD"], catalog["TELESCOPE"], catalog["FILE"])):
        if not mask[idx] or filename is None:  # no need to do anything if already bad
            mask[idx] = False
            continue
        result = visit_spectra(base_path, field, telescope, filename)
        # if file is missing locally, attempt to download
        if result == 1:
            mask[idx] = False
            continue
        result = combined_spectra(base_path, field, apogee_id, telescope)
        # if file is missing locally, attempt to download
        if result == 1:
            mask[idx] = False
    return mask


def download_allstar(base_path):
    fullfoldername = os.path.join(base_path, "/spectro/aspcap/dr17/synspec_rev1/")
    # Check if directory exists
    if not os.path.exists(fullfoldername):
        os.makedirs(fullfoldername)
    filename = "allStar-dr17-synspec_rev1.fits"
    fullfilename = os.path.join(fullfoldername, filename)
    url = f"https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/{filename}"
    urllib.request.urlretrieve(url, fullfilename)


def combined_spectra(
        base_path,
        field,
        apogee,
        telescope
        ):
    aspcap_code = "synspec_rev1"
    str1 = f"https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/{aspcap_code}/{telescope}/{field}/"

    filename = f"aspcapStar-dr17-{apogee}.fits"
    urlstr = str1 + filename

    # check folder existence
    fullfoldername = os.path.join(
        base_path,
        f"spectro/aspcap/dr17/{aspcap_code}/{telescope}",
        str(f"{field}"),
    )
    if not os.path.exists(fullfoldername):
        os.makedirs(fullfoldername)

    fullfilename = os.path.join(fullfoldername, filename)

    if os.path.exists(fullfilename):
        return fullfilename
    else:
        try:
            urllib.request.urlretrieve(urlstr, fullfilename)
            return 0
        except urllib.error.HTTPError as emsg:
            return 1  # error code


def visit_spectra(
    base_path,
    field,
    telescope,
    filename,
):
    str1 = f"https://data.sdss.org/sas/dr17/apogee/spectro/redux/dr17/stars/{telescope}/{field}/"
    urlstr = str1 + filename

    fullfoldername = os.path.join(
        base_path,
        f"spectro/redux/dr17/stars/{telescope}/",
        str(f"{field}"),
    )

    if not os.path.exists(fullfoldername):
        os.makedirs(fullfoldername)

    fullfilename = os.path.join(fullfoldername, filename)
    if os.path.exists(fullfilename):
        return fullfilename
    else:
        try:
            urllib.request.urlretrieve(urlstr, fullfilename)
            return 0
        except urllib.error.HTTPError as emsg:
            return 1  # error code


def processing_fn(raw_filename, continuum_filename):
    """Parallel processing function reading all requested spectra from one plate."""

    # Load the visit spectra file
    hdus = fits.open(raw_filename)
    raw_flux = hdus[1].data[0]
    raw_ivar = 1 / hdus[2].data[0] ** 2
    mask_spec = hdus[2].data[0]

    # Load the combined spectra file
    hdus = fits.open(continuum_filename)
    continuum_flux = hdus[1].data
    continuum_ivar = 1 / hdus[2].data ** 2

    # very rough estimate
    # https://www.sdss4.org/dr17/irspec/spectra/
    lsf_sigma = np.ones_like(raw_flux)
    lsf_sigma[:blue_end] *= 0.326
    lsf_sigma[blue_end:green_end] *= 0.283
    lsf_sigma[green_end:] *= 0.236

    # Return the results
    return {
        "spectrum_lambda": lam_cropped,
        "spectrum_flux": raw_flux,
        "spectrum_ivar": raw_ivar,
        # pixel level bitmask
        # see https://www.sdss4.org/dr17/irspec/apogee-bitmasks/#APOGEE_PIXMASK:APOGEEbitmaskforindividualpixelsinaspectrum
        'spectrum_lsf_sigma': lsf_sigma,
        "spectrum_bitmask": mask_spec,
        "pseudo_continuum_spectrum_flux": continuum_flux,
        "pseudo_continuum_spectrum_ivar": continuum_ivar,
    }


def save_in_standard_format(args):
    """
    This function takes care of iterating through the different input files
    corresponding to this healpix index, and exporting the data in standard format.
    """
    catalog, output_filename, apogee_data_path = args
    # Create the output directory if it does not exist
    if not os.path.exists(os.path.dirname(output_filename)):
        os.makedirs(os.path.dirname(output_filename))

    # Rename columns to match the standard format
    catalog["ra"] = catalog["RA"]
    catalog["dec"] = catalog["DEC"]
    catalog["object_id"] = catalog["APOGEE_ID"]
    catalog["radial_velocity"] = catalog["VHELIO_AVG"]
    catalog["restframe"] = np.ones(len(catalog), dtype=bool)

    # Preparing the arguments for the parallel processing
    # Process all files
    results = []
    for i in catalog:
        telescope = i["TELESCOPE"]
        field = i["FIELD"]
        filename = i["FILE"]
        apogee_id = i["APOGEE_ID"]
        results.append(
            processing_fn(
                visit_spectra(apogee_data_path, field, telescope, filename),
                combined_spectra(apogee_data_path, field, apogee_id, telescope)      
            )
        )

    # Aggregate all spectra into an astropy table
    spectra = Table({k: np.vstack([d[k] for d in results]) for k in results[0].keys()})

    # crop apogee spectra
    spectra["spectrum_flux"] = spectra["spectrum_flux"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]
    spectra["spectrum_ivar"] = spectra["spectrum_ivar"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]
    spectra["spectrum_bitmask"] = spectra["spectrum_bitmask"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]    
    spectra["spectrum_lsf_sigma"] = spectra["spectrum_lsf_sigma"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]
    spectra["pseudo_continuum_spectrum_flux"] = spectra["pseudo_continuum_spectrum_flux"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]    
    spectra["pseudo_continuum_spectrum_ivar"] = spectra["pseudo_continuum_spectrum_ivar"][
        :, np.r_[blue_start:blue_end, green_start:green_end, red_start:red_end]
    ]
    # Join on target id with the input catalog
    # catalog = join(catalog, spectra, keys="object_id", join_type="inner")
    catalog = hstack([catalog, spectra])

    # Making sure we didn't lose anyone
    assert (
        len(catalog) == len(spectra)
    ), "There was an error in the join operation, probably some spectra files are missing"

    # Save all columns to disk in HDF5 format
    with h5py.File(output_filename, "w") as hdf5_file:
        for key in catalog.colnames:
            hdf5_file.create_dataset(key.lower(), data=catalog[key])
    return 1


def main(args):
    # Load the catalog file and apply main cuts
    path_to_read = os.path.join(
            args.apogee_data_path,
            "spectro/aspcap/dr17/synspec_rev1",
            "allStar-dr17-synspec_rev1.fits",
        )
    if not os.path.exists:
        download_allstar(args.apogee_data_path)
    catalog = Table.read(path_to_read, hdu=1)

    catalog = catalog[selection_fn(args.apogee_data_path, catalog)]

    # Add healpix index to the catalog
    catalog["healpix"] = hp.ang2pix(
        64, catalog["RA"].filled(), catalog["DEC"].filled(), lonlat=True, nest=True
    )
    catalog = catalog.group_by(["healpix"])

    # Preparing the arguments for the parallel processing
    map_args = []
    for group in catalog.groups:
        # Create a filename for the group
        group_filename = os.path.join(
            args.output_dir,
            "apogee/healpix={}/001-of-001.hdf5".format(group["healpix"][0]),
        )
        map_args.append((group, group_filename, args.apogee_data_path))

    # Run the parallel processing
    with Pool(args.num_processes) as pool:
        results = list(
            tqdm(pool.imap(save_in_standard_format, map_args), total=len(map_args))
        )

    if sum(results) != len(map_args):
        print(
            "There was an error in the parallel processing, some files may not have been processed correctly"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download (if not existed already) and Extracts spectra from all SDSS spectra"
    )
    parser.add_argument(
        "--apogee_data_path",
        type=str,
        help="Path to the local copy of the APOGEE data where path should end with path /dr17/apogee/",
    )
    parser.add_argument("--output_dir", type=str, help="Path to the output directory")
    parser.add_argument(
        "--num_processes",
        type=int,
        default=10,
        help="The number of processes to use for parallel processing",
    )
    args = parser.parse_args()

    main(args)
