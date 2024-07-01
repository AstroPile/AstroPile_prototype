import argparse
from functools import partial
import os
import urllib.request

from tqdm.contrib.concurrent import process_map


def _download_file(f, args):
    f = f.strip()
    filename = f.split("/")[-1]
    os.makedirs(args.output_dir, exist_ok=True)
    savename = os.path.join(args.output_dir, filename)
    if not os.path.exists(savename):
        urllib.request.urlretrieve(f, savename)


def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    with open("file_list.txt") as f:
        files = f.readlines()

    if args.tiny:
        files = files[:1]

    if args.aria2:
        files = " ".join([f'"{f.strip()}"' for f in files])
        os.system(f"aria2c -j16 -s16 -x16 -c -d {args.output_dir} -Z {files}")

    else:
        _download = partial(_download_file, args=args)

        process_map(_download, files, max_workers=16, chunksize=1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tiny",
        help="download a single file only",
        action="store_true",
    )
    parser.add_argument(
        "--aria2",
        help="use aria2c to download files",
        action="store_true",
    )
    parser.add_argument(
        "--output_dir",
        help="output directory",
        default=".",
    )
    args = parser.parse_args()
    main(args)
