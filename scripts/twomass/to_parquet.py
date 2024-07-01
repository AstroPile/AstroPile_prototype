import argparse
from functools import partial
import glob
import os

import healpy as hp
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq
from tqdm.contrib.concurrent import process_map

from twomass import _mapping


def ang2pix(nside, ra, dec):
    return hp.ang2pix(nside, ra, dec, lonlat=True, nest=True)


def read_table(filename, args):
    table = pcsv.read_csv(
        filename,
        read_options=pcsv.ReadOptions(
            use_threads=True, column_names=list(_mapping.keys())
        ),
        parse_options=pcsv.ParseOptions(delimiter="|"),
        convert_options=pcsv.ConvertOptions(column_types=_mapping, null_values=["\\N"]),
    )
    healpix = ang2pix(args.nside, table["ra"], table["decl"])
    table = table.append_column("healpix", [healpix])
    return table


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--nside", type=int, default=16)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--file_prefix", type=str, required=True)
    parser.add_argument("--tiny", action="store_true", default=False)
    args = parser.parse_args()

    filenames = glob.glob(f"{args.data_dir}/{args.file_prefix}*")

    results = process_map(partial(read_table, args=args), filenames)
    table = pa.concat_tables(results)  # , promote_options="permissive")
    if args.tiny:
        table = table.take(list(range(2500)))

    os.makedirs(args.output_dir, exist_ok=True)
    pq.write_to_dataset(table, args.output_dir, partition_cols=["healpix"])
