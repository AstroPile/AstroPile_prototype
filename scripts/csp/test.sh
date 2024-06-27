#!/bin/bash

# Download the dataset
if python download_data.py ./; then
    echo "Download CSP dataset successful"
else
    echo "Download CSP dataset failed"
    exit 1
fi

# First build the parent sample and save both raw and H5 to current directory
if python build_parent_sample.py ./CSPDR3 ./ --tiny; then
    echo "Build parent sample for CSP successful"
else
    echo "Build parent sample for CSP failed"
    exit 1
fi

# Try to load the dataset with hugging face dataset
if python -c "from datasets import load_dataset; dset = load_dataset('./csp.py', trust_remote_code=True, split='train'); print(f'loaded dataset with {len(dset)} examples'); next(iter(dset));"; then
    echo "Load dataset for CSP successful"
else
    echo "Load dataset for CSP failed"
    exit 1
fi
