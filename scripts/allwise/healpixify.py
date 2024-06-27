import healpy as hp
import pyarrow.parquet as pq


def ang2pix(ra, dec):
    return hp.ang2pix(nside=args.nside, theta=ra, phi=dec, lonlat=True, nest=True)


def main(args):
    table = pq.read_table(args.input_dir)
    hpix = ang2pix(table["ra"], table["dec"])
    table = table.drop_columns(["healpix_k0", "healpix_k5"])
    table = table.append_column("healpix", [hpix])
    pq.write_to_dataset(table, args.output_dir, partition_cols=["healpix"])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir",
        type=str,
        help="Path to the local data directory",
    )
    parser.add_argument("--output_dir", type=str, help="Path to the output directory")
    parser.add_argument("--nside", type=int, help="nside for healpix")
    args = parser.parse_args()

    main(args)
