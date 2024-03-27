import os

import h5py
import healpy as hp
import numpy as np


def save_in_standard_format(args):
    source_file, output_filename, selection_mask, healpix_id = args
    catalog = h5py.File(source_file, "r")

    if not os.path.exists(os.path.dirname(output_filename)):
        os.makedirs(os.path.dirname(output_filename))

    with h5py.File(output_filename, "w") as f:
        # Save the data
        for key in catalog.keys():
            f.create_dataset(key, data=catalog[key][selection_mask])
        f.create_dataset(
            "healpix", data=np.repeat(healpix_id, len(selection_mask)).astype(np.int64)
        )

    catalog.close()

    return 1


def ang2pix(ra, dec):
    return hp.ang2pix(nside=args.nside, theta=ra, phi=dec, lonlat=True, nest=True)


def main(args):
    from tqdm.auto import tqdm
    from multiprocessing import Pool

    source_file = args.input_file
    catalog = h5py.File(source_file, "r")

    with Pool(args.num_procs) as pool:
        healpix = np.concatenate(
            list(
                tqdm(
                    pool.starmap(
                        ang2pix,
                        zip(
                            np.array_split(np.array(catalog["ra"]), args.num_procs),
                            np.array_split(np.array(catalog["dec"]), args.num_procs),
                        ),
                    ),
                    total=args.num_procs,
                )
            )
        )

    hp_groups = np.unique(healpix)

    catalog.close()

    # Preparing the arguments for the parallel processing
    map_args = []
    for hp_ix in hp_groups:
        # Create a filename for the group
        output_filename = os.path.join(
            args.output_dir, f"gaia/healpix={hp_ix}/001-of-001.hdf5"
        )
        selection_mask = np.where(healpix == hp_ix)[0]
        map_args.append((source_file, output_filename, selection_mask, hp_ix))

    # Run the parallel processing
    with Pool(args.num_procs) as pool:
        results = list(
            tqdm(pool.imap(save_in_standard_format, map_args), total=len(map_args))
        )

    if sum(results) != len(map_args):
        print(
            "There was an error in the parallel processing, some files may not have been processed correctly"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file",
        type=str,
        help="Path to the local data file",
    )
    parser.add_argument("--output_dir", type=str, help="Path to the output directory")
    parser.add_argument("--nside", type=int, help="nside for healpix")
    parser.add_argument(
        "--num_procs",
        type=int,
        default=10,
        help="The number of processes to use for parallel processing",
    )
    args = parser.parse_args()

    main(args)
