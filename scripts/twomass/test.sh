#!/bin/bash

echo ==========
echo 2MASS TEST
echo ==========

echo ========================
echo DOWNLOADING TINY DATASET
echo ========================

# download files
python3 download_parts.py --tiny --output_dir _2mass/psc --aria2

echo =====================================
echo PREPARING PARTITIONED PARQUET DATASET
echo =====================================

# healpixify 
python3 to_parquet.py --nside 16 --data_dir _2mass/psc --output_dir _2mass_pq/psc --file_prefix psc --tiny

echo ==================
echo CONVERTING TO HDF5 
echo ==================

# convert
python3 to_hdf5.py --data_dir _2mass_pq --output_dir 2mass

echo ==================
echo TESTING HF LOADING
echo ==================

# test loading
python3 test_load.py

echo ============
echo TEST PASSED! 
echo ============
