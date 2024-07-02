import multiprocessing
import h5py
import argparse
from glob import glob
from tqdm import tqdm


def check_file(filename: str):
    try:
        with h5py.File(filename, "r") as f:
            f["object_id"][:]
    except Exception as e:
        return filename, e


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tests whether all hdf5 can be read correctly for a given sample."
    )
    parser.add_argument("sample_path", type=str, help="Path of the sample")
    parser.add_argument("--n_proc", type=int, default=1, help="Number of processes")
    parser.add_argument("--chunks", type=int, default=1, help="Chunk size per process")
    args = parser.parse_args()
    root = args.sample_path
    n_proc = args.n_proc
    chunks = args.chunks
    path = f"{root}/*/*/*.hdf5"
    files = glob(path)
    n_problematic_files = 0
    with multiprocessing.Pool(n_proc) as pool:
        results = tqdm(
            pool.imap_unordered(check_file, files, chunksize=chunks), total=len(files)
        )
        for r in results:
            if r is not None:
                filename, error = r
                print(f"Fail to load {filename} due to error {error}")
                n_problematic_files += 1
    print(f"Failed to load {n_problematic_files} files out of {len(files)}")
