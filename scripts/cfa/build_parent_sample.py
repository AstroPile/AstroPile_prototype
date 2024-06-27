import argparse
import os
import shutil

import h5py
import healpy as hp
import numpy as np
import pandas as pd


def get_str_dtype(arr):
    str_max_len = int(np.char.str_len(arr).max())
    return h5py.string_dtype(encoding="utf-8", length=str_max_len)


def convert_dtype(arr):
    if np.issubdtype(arr.dtype, np.floating):
        dtype = np.float32
    elif np.issubdtype(arr.dtype, np.str_):
        dtype = get_str_dtype(arr)
    else:
        dtype = arr.dtype
    return arr.astype(dtype)


def cfa_snII_bpf(file_dir, data, metadata, keys_data, keys_metadata, tiny=False, **kwargs):
    info = {}
    with open("CFA_SNII_COORDS.txt", "r") as f:
        for line in f.readlines():
            info[line.split()[0]] = line.split()[1:]
    optical_df = pd.read_csv(
        os.path.join(file_dir, "STDSYSTEM_LC.txt"),
        comment="#",
        sep=r"\s+",
        names=["name", "FLT", "time", "N", "mag", "mag_err"],
    )
    nir_df = pd.read_csv(
        os.path.join(file_dir, "NIR_LC.txt"),
        comment="#",
        sep=r"\s+",
        names=["name", "FLT", "time", "mag", "mag_err"],
    )
    df = pd.concat([optical_df, nir_df])
    unique_names = set(df["name"])
    if tiny:
        unique_names = list(unique_names)[:10]
    for sn_name in unique_names:
        mask = np.where(df["name"] == sn_name)
        for key in keys_data:
            data[key].append(np.array(df[key])[mask])
        metadata["object_id"].append(sn_name)
        metadata["ra"].append(float(info[sn_name][0]))
        metadata["dec"].append(float(info[sn_name][1]))
        metadata["redshift"].append(0)
        # Assuming all are SNe II. There may be unlabeled subtypes.
        metadata["obj_type"].append(" ".join(info[sn_name][2:]))

    num_examples = len(metadata["object_id"])

    return num_examples, data, metadata


def cfa3_bpf(file_dir, data, metadata, keys_data, keys_metadata, tiny=False, **kwargs):
    bandpass_dict = {
        "1": "U",
        "2": "B",
        "3": "V",
        "4": "R",
        "5": "I",
        "13": "r'",
        "14": "i'",
    }
    info = {}
    with open("CFA3_COORDS.txt", "r") as f:
        for line in f.readlines():
            info[line.split()[0]] = line.split()[1:]
    f = open(os.path.join(file_dir, "cfa3lightcurves.standardsystem.txt"), "r")
    current_sn = ""
    data_ = dict(zip(keys_data, ([] for _ in keys_data)))
    for line in f.readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sn"):
            # Do not add data_ on first sn
            if current_sn != "":
                for key in keys_data:
                    data[key].append(np.array(data_[key]))
            current_sn = line.strip("\n")
            data_ = dict(zip(keys_data, ([] for _ in keys_data)))
            sn_name = "SN20" + current_sn[2:]
            metadata["object_id"].append(sn_name)
            metadata["obj_type"].append(" ".join(info[sn_name][2:]))
            metadata["ra"].append(float(info[sn_name][0]))
            metadata["dec"].append(float(info[sn_name][1]))
            metadata["redshift"].append(0)
        else:
            bp, mjd, mag, dmag = line.split()
            data_["FLT"].append(bandpass_dict[bp])
            data_["time"].append(float(mjd))
            data_["mag"].append(float(mag))
            data_["mag_err"].append(float(dmag))
    # last sn ends with a data line. Add data_ manually
    for key in keys_data:
        data[key].append(np.array(data_[key]))
    f.close()
    num_examples = len(metadata["object_id"])
    return num_examples, data, metadata


def cfa_generic_bpf(file_dir, data, metadata, keys_data, keys_metadata, dataset, tiny=False, **kwargs):
    if dataset == "cfa_SECCSN":
        info_file = "CFA_SECCSN_COORDS.txt"
        file_name = "lc.standardsystem.sesn_allphot.dat"
        columns = ["name", "FLT", "time", "mag", "mag_err", "survey"]
    elif dataset == "cfa4":
        info_file = "CFA4_COORDS.txt"
        file_name = "cfa4.lc.stdsystem.fi.ascii"
        columns = [
            "name",
            "FLT",
            "time",
            "N",
            "sigma_pipe",
            "sigma_phot",
            "mag",
            "mag_err",
        ]
    info = {}
    with open(info_file, "r") as f:
        for line in f.readlines():
            info[line.split()[0]] = line.split()[1:]
    df = pd.read_csv(
        os.path.join(file_dir, file_name),
        comment="#",
        sep=r"\s+",
        names=columns,
    )
    unique_names = set(df["name"])
    if tiny:
        unique_names = list(unique_names)[:10]
    for sn_name in unique_names:
        mask = np.where(df["name"] == sn_name)
        for key in keys_data:
            data[key].append(np.array(df[key])[mask])
        metadata["object_id"].append("SN" + sn_name)
        if not sn_name.startswith("snf"):
            metadata["ra"].append(float(info["SN" + sn_name][0]))
            metadata["dec"].append(float(info["SN" + sn_name][1]))
        else:
            metadata["ra"].append(0)
            metadata["dec"].append(0)
        metadata["obj_type"].append(" ".join(info["SN" + sn_name][2:]))
        metadata["redshift"].append(0)

    num_examples = len(metadata["object_id"])

    return num_examples, data, metadata


survey_specific_logic = {
    "cfa3": cfa3_bpf,
    "cfa_SECCSN": cfa_generic_bpf,
    "cfa4": cfa_generic_bpf,
    "cfa_snII": cfa_snII_bpf,
}


def main(args):
    # Retrieve file paths
    file_dir = args.data_path
    output_dir = args.output_dir
    dataset = args.dataset
    tiny = args.tiny
    dirty = args.dirty

    keys_data = ["time", "mag", "mag_err", "FLT"]
    keys_metadata = ["object_id", "redshift", "ra", "dec", "obj_type"]
    data = dict(zip(keys_data, ([] for _ in keys_data)))
    metadata = dict(zip(keys_metadata, ([] for _ in keys_metadata)))
    num_examples, data, metadata = survey_specific_logic[dataset](
        file_dir,
        data,
        metadata,
        keys_data,
        keys_metadata,
        dataset=dataset,
        tiny=tiny,
    )

    # Create an array of all bands in the dataset
    all_bands = np.unique(np.concatenate(data["FLT"]))

    # Remove band from keys_data as the timeseries will be arranged by band
    keys_data.remove("FLT")

    for i in range(num_examples):
        # For this example, find the band with the most observations
        # and store the number of observations as max_length
        _, count = np.unique(data["FLT"][i], return_counts=True)
        max_length = count.max()

        # Create mask to select data from each timeseries by band
        mask = np.expand_dims(all_bands, 1) == data["FLT"][i]

        for key in keys_data:
            timeseries_all_bands = (
                []
            )  # Stores a particular timeseries (as specified by the key) in all bands
            for j in range(len(all_bands)):
                timeseries_band = data[key][i][
                    mask[j]
                ]  # Select samples from timeseries for a particular band
                timeseries_band = np.pad(  # Pad single band timeseries to max_length
                    timeseries_band,
                    (0, max_length - len(timeseries_band)),
                    mode="constant",
                    constant_values=-99 if key == "MJD" else 0,
                )
                timeseries_all_bands.append(timeseries_band)
            timeseries_all_bands = convert_dtype(np.array(timeseries_all_bands))
            data[key][i] = timeseries_all_bands

    # Convert metadata to numpy arrays and cast to required datatypes
    for key in keys_metadata:
        metadata[key] = convert_dtype(np.array(metadata[key]))

    # Add healpix to metadata
    keys_metadata.append("healpix")
    metadata["healpix"] = hp.ang2pix(
        16, metadata["ra"], metadata["dec"], lonlat=True, nest=True
    )

    # Cast bands to required datatype
    all_bands = convert_dtype(all_bands)

    # Establish conversions to standard names
    keys_all = keys_metadata + keys_data
    name_conversion = dict(zip(keys_all, keys_all))

    # Make output directories labelled by healpix
    unique_healpix = np.unique(metadata["healpix"])
    healpix_num_digits = len(str(hp.nside2npix(16)))
    for healpix in unique_healpix:
        healpix = str(healpix).zfill(healpix_num_digits)
        os.makedirs(os.path.join(output_dir, f"healpix={healpix}"), exist_ok=True)

    # Save data as hdf5 grouped into directories by healpix
    for i in range(num_examples):
        healpix = str(metadata["healpix"][i]).zfill(healpix_num_digits)
        object_id = metadata["object_id"][i].decode("utf-8")
        path = os.path.join(
            output_dir, f"healpix={healpix}", f"example_{object_id}.hdf5"
        )
        with h5py.File(path, "w") as hdf5_file:
            # Save metadata
            for key in keys_metadata:
                hdf5_file.create_dataset(name_conversion[key], data=metadata[key][i])
            # Save bands
            hdf5_file.create_dataset("bands", data=all_bands)
            # Save timeseries
            for key in keys_data:
                hdf5_file.create_dataset(name_conversion[key], data=data[key][i])

    # Remove original data (data has now been reformatted and saved as hdf5)
    if not dirty:
        shutil.rmtree(file_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract data and convert to standard time-series data format."
    )
    parser.add_argument(
        "data_path", type=str, help="Path to the local copy of the data"
    )
    parser.add_argument("output_dir", type=str, help="Path to the output directory")
    parser.add_argument("dataset", type=str, help="Dataset to be prepared")
    parser.add_argument(
        "--tiny", action="store_true", help="Use a small subset of the data for testing"
    )
    parser.add_argument(
        "--dirty", action="store_true", help="Do not remove the original data"
    )
    args = parser.parse_args()

    main(args)
