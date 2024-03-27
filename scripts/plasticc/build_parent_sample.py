from pathlib import Path
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from multiprocessing import Pool
from tqdm import tqdm
import h5py
import pandas as pd
import pdb

def save_in_standard_format(args):
    """ This function iterates through an input metadata/lightcurve data pair and saves the data in a standard format.
    """
    metadata, lcdata, output_filename = args
    if output_filename.exists():
        print(f"reading existing data from {output_filename}")
        return 1

    metadata = metadata[metadata["object_id"].isin(lcdata["object_id"].unique())]
    print(len(metadata.object_id.unique()), len(lcdata.object_id.unique()), flush=True)
    print(len(metadata.object_id.unique()), len(lcdata.object_id.unique()))

    # Create the output directory if it does not exist
    if not Path(output_filename).parent.exists():
        Path(output_filename).parent.mkdir(parents=True)

    # find longest light curve
    max_length = lcdata.groupby('object_id').size().max()

    # process this file
    ids = np.unique(lcdata['object_id'])
    objects = []
    for i, object_id in enumerate(ids):
        if i % 50_000 == 0:
            print(f"{output_filename.name}: {i}/{len(ids)}", flush=True)
        obj_metadata = metadata[metadata['object_id'] == object_id]
        lcdata = lcdata[lcdata['object_id'] == object_id]

        lc = np.zeros((6, 3, max_length))
        for band in range(5):
            band_data = lcdata[lcdata['passband'] == band]
            lc[band, 0, :len(band_data)] = band_data['mjd']
            lc[band, 1, :len(band_data)] = band_data['flux']
            lc[band, 2, :len(band_data)] = band_data['flux_err']

        objects.append({
            'object_id': object_id,
            'ra': obj_metadata['ra'].values[0],
            'dec': obj_metadata['decl'].values[0],
            'hostgal_specz': obj_metadata['hostgal_specz'].values[0],
            'hostgal_photoz': obj_metadata['hostgal_photoz'].values[0],
            'redshift': obj_metadata['true_z'].values[0],
            'obj_type': obj_metadata['true_target'].values[0],
            'lightcurve': lc,
        })

    print(f"{output_filename.name}: finished processing", flush=True)
    objects = Table({k: [o[k] for o in objects] for k in objects[0].keys()})
    print(f"{output_filename.name}: finished creating table", flush=True)
    assert len(objects) == len(ids), f"Lost objects during preprocessing: {len(objects)} != {len(ids)}"

    # Save all columns to disk in HDF5 format
    with h5py.File(output_filename, 'w') as hdf5_file:
        for key in objects.colnames:
            hdf5_file.create_dataset(key, data=objects[key])
    print(f"{output_filename.name}: finished writing h5", flush=True)
    return 1

def download_plasticc_data(output_path):
    from urllib.request import urlretrieve

    if not Path(output_path).exists():
        Path(output_path).mkdir(parents=True)

    url_base = "https://zenodo.org/records/2539456/files/{}?download=1"
    urls = [
            url_base.format("plasticc_train_metadata.csv.gz"),
            url_base.format("plasticc_train_lightcurves.csv.gz"),
            url_base.format("plasticc_test_metadata.csv.gz"),
            url_base.format("plasticc_test_lightcurves_01.csv.gz"),
            url_base.format("plasticc_test_lightcurves_02.csv.gz"),
            url_base.format("plasticc_test_lightcurves_03.csv.gz"),
            url_base.format("plasticc_test_lightcurves_04.csv.gz"),
            url_base.format("plasticc_test_lightcurves_05.csv.gz"),
            url_base.format("plasticc_test_lightcurves_06.csv.gz"),
            url_base.format("plasticc_test_lightcurves_07.csv.gz"),
            url_base.format("plasticc_test_lightcurves_08.csv.gz"),
            url_base.format("plasticc_test_lightcurves_09.csv.gz"),
            url_base.format("plasticc_test_lightcurves_10.csv.gz"),
            url_base.format("plasticc_test_lightcurves_11.csv.gz"),
        ]
    for url in tqdm(urls):
        filename = url.split('/')[-1].split('?')[0]
        filepath = Path(output_path) / filename
        if filepath.exists():
            print(f"reading existing plasticc data from {filename}")
        else:
            print(f"no {filename} found, downloading ...")
            urlretrieve(url, (Path(output_path) / filename).resolve())
    print("done!")

def main(args):
    # Load PLAsTiCC data locally if exists, else download from Zenodo
    download_plasticc_data(args.plasticc_data_path)

    print("Rewriting training data into standard format...")
    # process training data
    metadata = pd.read_csv(Path(args.plasticc_data_path) / "plasticc_train_metadata.csv.gz")
    lcdata = pd.read_csv(Path(args.plasticc_data_path) / "plasticc_train_lightcurves.csv.gz")
    save_in_standard_format((metadata, lcdata, Path(args.output_path) / "plasticc_train.hdf5"))

    # process test data
    print("Rewriting test data into standard format...")
    metadata = pd.read_csv(Path(args.plasticc_data_path) / "plasticc_test_metadata.csv.gz")
    map_args = [[metadata, pd.read_csv(Path(args.plasticc_data_path) / f"plasticc_test_lightcurves_{i:02d}.csv.gz"), Path(args.output_path) / f"plasticc_test_{i:02d}.hdf5"] for i in range(1, 12)]

    # Run the parallel processing
    with Pool(args.num_processes) as pool:
        results = list(tqdm(pool.imap(save_in_standard_format, map_args), total=len(map_args)))

    if sum(results) != len(map_args):
        print("There was an error in the parallel processing, some files may not have been processed correctly")

    # clean up the original data files
    print("Cleaning up original data files...")
    for i in range(1, 12):
        Path(args.plasticc_data_path / f"plasticc_test_lightcurves_{i:02d}.csv.gz").unlink()
    Path(args.plasticc_data_path / "plasticc_train_metadata.csv.gz").unlink()
    Path(args.plasticc_data_path / "plasticc_train_lightcurves.csv.gz").unlink()
    Path(args.plasticc_data_path / "plasticc_test_metadata.csv.gz").unlink()

    print("All done!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts lightcurves and metadata from all PLAsTiCC data')
    parser.add_argument('plasticc_data_path', type=str, help='Path to the local copy of the PLAsTiCC data')
    parser.add_argument('output_path', type=str, help='Path to the output directory')
    parser.add_argument('--num_processes', type=int, default=10, help='The number of processes to use for parallel processing')
    args = parser.parse_args()

    main(args)
