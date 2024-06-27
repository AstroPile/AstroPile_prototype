#!/bin/bash

echo =============
echo AllWISE TEST
echo =============

echo ========================
echo DOWNLOADING TINY DATASET
echo ========================

# download files
python3 download_parts.py --tiny --output_dir _allwise

echo ===================
echo HEALPIXIFYING FILES
echo ===================

# healpixify 
python3 healpixify.py --nside 16 --input_dir _allwise --output_dir _allwise_hp

echo ==================
echo CONVERTING TO HDF5 
echo ==================

# convert
python3 to_hdf5.py --data_dir _allwise_hp --output_dir allwise

echo ==================
echo TESTING HF LOADING
echo ==================

# test loading
python3 test_load.py

echo ============
echo TEST PASSED! 
echo ============
