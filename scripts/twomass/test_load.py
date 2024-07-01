import datasets

ds_hdf5 = datasets.load_dataset(
    "twomass.py",
    trust_remote_code=True,
    split="train",
    streaming=True,
)

batch_hdf5 = next(iter(ds_hdf5))

ds_parquet = datasets.load_dataset(
    "_2mass_pq",
    trust_remote_code=True,
    split="train",
    streaming=True,
)

batch_parquet = next(iter(ds_parquet))
