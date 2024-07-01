# 2MASS

This folder contains a dataset for the Two Micron All Sky Survey (2MASS), which provides all-sky near-IR photometry in the J, H, and Ks bands at 1.235, 1.662, and 2.159 μm. Current there is support for the Point Source Catalog, which consists of over 500 million stars and galaxies.

For more information on 2MASS, see https://www.ipac.caltech.edu/2mass/overview/about2mass.html and https://www.ipac.caltech.edu/2mass/releases/allsky/doc/explsup.html.

For more information on the Point Source Catalog, see https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2.html, and particular https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2a.html for details on the format and columns.

## Data Preparation

The data can be downloaded using the `download_parts.py` script, which will grab the 2MASS PSC data in gzipped csv format. You can then run the `to_parquet.py` script, which will read in all the gzipped csv files, turn them into PyArrow tables, add a healpix column for partitioning, and then merge everything into a partitioned Parquet dataset. This can be used with HF dataset directly. Running `to_hdf5.py`, will convert the Parquet files into HDF5 files, after which you can use the `twomass.py` descriptor to load the files into HF datasets. As an example of the full preparation, you can see the `test.sh` script.

## Dataset

See https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2a.html for details on the format and columns. All columns are exposed.
