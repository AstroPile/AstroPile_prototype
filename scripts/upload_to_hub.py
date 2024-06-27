# Small script to run the dataset preparation for all astropile datasets
from datasets import load_dataset, Dataset
from huggingface_hub import DatasetCardData, DatasetCard
from glob import glob
from functools import partial
import pickle
import sys
import shutil

# Define the destination folder where the scripts will be copied
destination_folder = "/mnt/ceph/users/polymathic/MultimodalUniverse"

# Max size of the dataset to upload
max_size = 100 # GB

only_cards = False

# Listing all the datasets in the target directorty
datasets = glob('/mnt/ceph/users/polymathic/MultimodalUniverse/*/*.py')
datasets = [d.split('/')[-1].split('.py')[0] for d in datasets]

print("Found the following datasets:")
print(datasets)

def gen_from_iterable_dataset(iterable_ds):
    yield from iterable_ds

for dataset in datasets:
    print(f'Preparing {dataset} dataset') 

    # Copy the file to the destination folder
    shutil.copyfile(f"./{dataset}/{dataset}.py",
                    f"{destination_folder}/{dataset}/{dataset}.py")
    print("Copied file:", f"{dataset}/{dataset}.py")

    try:
        dset = load_dataset(f'{destination_folder}/{dataset}', 
                            trust_remote_code=True, streaming=True, split='train')

        # Extract the first element of the dataset to evaluate the size
        top = next(iter(dset))
        element_size = sys.getsizeof(pickle.dumps(top))
        max_elements = int(max_size * 1e9 / element_size)
        # TODO: Remove this line and instead use the max size to know how many elements to upload
        max_elements = 1000
        print(f"Element size: {element_size} bytes, so we can upload {max_elements} elements to the hub.")

        if not only_cards:
            # Restrict the dataset to the maximum number of elements as to remain within budget
            dset = dset.take(max_elements)

            ds = Dataset.from_generator(partial(gen_from_iterable_dataset, dset), 
                                        features=dset.features,
                                        cache_dir=f'/tmp/hf_cache/{dataset}')

            # Save the dataset
            ds.push_to_hub(f'MultimodalUniverse/{dataset}', config_name=dset.config_name)
            print(f"Dataset {dataset} uploaded to the hub")

        # Create the dataset card
        card = DatasetCardData(
            description=dset.description,
            homepage=dset.homepage,
            version=dset.version.__str__(),
            citation=dset.citation,
        )
        content=f"""
---
{card.to_yaml()}
---

# {dataset.capitalize()} Dataset

{dset.license}

{dset.description}

{dset.citation}
"""
        c = DatasetCard(content)
        c.push_to_hub(f'MultimodalUniverse/{dataset}')

    except Exception as e:
        print(f"Failed to upload {dataset} dataset")
        print(e)
        continue
