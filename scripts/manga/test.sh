#!/bin/bash

# install the sdss-access tool
pip install sdss_access

# download a small sample
echo "Downloading sample MaNGA data"
if python access_transfer.py --limit --destination_path .; then
    echo "Downloading data succeeded"
    cp dr17/manga/spectro/analysis/v3_1_1/3.1.0/dapall*.fits .
    cp dr17/manga/spectro/redux/v3_1_1/drpall*.fits .
else
    echo "Downloading data failed"
    exit 1
fi


# build the parent sample
echo "Building parent hdf5 sample"
if python build_parent_sample.py --manga_data_path . --output_dir . --num_processes 1 --tiny; then
    echo "Building parent sample for manga succeeded"
else
    echo "Building parent sample for manga failed"
    exit 1
fi

# try to load the dataset
echo "Loading dataset"
if python -c "from datasets import load_dataset; dset = load_dataset('./manga.py', trust_remote_code=True, split='train', streaming=True).with_format('numpy'); next(iter(dset));"; then
    echo "Loaded manga dataset succeeded"
else
    echo "Loading manga dataset failed"
    exit 1
fi