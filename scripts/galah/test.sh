#!/bin/bash

# First build the parent sample and save both raw and H5 to current directory
if python build_parent_sample.py "./downloaded_galah" "." --num_processes 1 --tiny; then
    echo "Build parent sample for GALAH successful"
else
    echo "Build parent sample for GALAH failed"
    exit 1
fi

# Try to load the dataset with hugging face dataset
if python -c "from datasets import load_dataset; dset = load_dataset('./galah.py', trust_remote_code=True, split='train').with_format('numpy'); print(f'loaded dataset with {len(dset)} examples'); next(iter(dset));"; then
    echo "Load dataset for GALAH successful"
else
    echo "Load dataset for GALAH failed"
    exit 1
fi
