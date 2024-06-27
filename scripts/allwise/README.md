# AllWISE

This folder contains a dataset for the AllWISE Source Catalog, which provides all-sky mid-IR photometry in the WISE W1, W2, W3, and W4 bands at 3.4, 4.6, 12 and 22 μm, along with position and apparent motion information for 747 million sources.

For more information on AllWISE, see https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec1_1.html.

For more information on the Source Catalog, see https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1.html, and particular https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html for details on the format and columns.

## Data Preparation

The data can be downloaded using the `download_parts.py` script, which will grab the AllWISE data in Parquet format. This is already partitioned into a healpix k=5 (nside=32) Parquet dataset, which you can immediately load into Huggingface datasets with `datasets.load_dataset(...)`. However, for consistency with the rest of the datasets in AstroPile, we prepare it into HDF5 files partitioned by healpix at k=4 (nside=16). This is done by running `healpixify.py`, which repartitions the Parquet files, and then `to_hdf5.py`, which will convert the Parquet files into HDF5 files. As an example of the full preparation, you can see the `test.sh` script.

## Dataset

See https://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html for details on the format and columns. All columns are exposed.
