# Design Document for the AstroPile


## Goals and Scope

The AstroPile aims to build a large dataset of diverse astronomical data that may be used to build
astronomical foundation models.

The current types of data that are being collected are:
  - Multiband images from different surveys and instruments
  - Spectra from different instruments

Contrary to previous datasets collected for machine learning applications, the AstroPile does not
format all data modalities into a uniform way, but rather collects sufficient and necessary metadata
to allow machine learning models to understand the context of each observation.

The AstroPile also includes pairing between different data modalities for the same astronomical objects.


## Dataset Formats

### Images

The organization of the data is based on the following strategy:

  - Each source survey is packaged as a separare parent dataset.
  - An `Example` in these parent datasets consists of all observations of an astronomical object from that survey. An example for the format of an `Example` in a survey like HSC would be:
  ```
    {
     "provenance": {"survey": "HSC", release: "1.0"},
     'observation_id': 111111,
     'ra': 124.,
     'dec': 0.,
     'healpix': 1234,
     'image': [{'array': array(224, 224), 'psf_fwhm': 0.7, 'pixel_size': 0.168, 'noise_std':  0.1, 'filter': 'r', 'units': 'nanomaggies', extinction: 0.1},
               {'array': array(224, 224), 'psf_fwhm': 0.7, 'pixel_size': 0.168, 'noise_std':  0.1, 'filter': 'i', 'units': 'nanomaggies', extinction: 0.1}],
      ...
      'forced.g_cmodel_mag': 22.5,
      'forced.g_cmodel_mag_err': 0.1,
      'forced.i_cmodel_mag': 21.5,
      'forced.i_cmodel_mag_err': 0.1,
      ...
    }
  ```

### Spectra

Spectral data should be organized with the following fields:

```
{
  "provenance": {"project": "DESI", "survey": "DESI", release: "DR1"},
  "observation_id": 111111,
  "ra": 124.,
  "dec": 0.,
  'flux': array(5000),
  'ivar': array(5000),
  'lsf_sigma': array(5000),
  'lambda': array(5000),
  'z': 0.01,
  'z_err': 1e-5,
  'extra': {}
}
```

### IFU Data

IFU data should be organized as the following example:
```
{
  "provenance": {"project": "SDSS", "survey": "MaNGA", release: "DR17"},
  "observation_id": 111111,
  "ra": 124.,
  "dec": 0.,
  "z": 0.01,
  "z_err": 1e-5,
  "spaxel_size": 0.5",
  "spaxels": [{'flux': array(5000), 'ivar': array(5000), 'lsf_sigma': array(5000), 'lambda': array(5000), 'x': 0, 'y': 0, 'extra': {}},
  ...
  ],
  "images": [
    {"label": "reconstructed_u", "filter": "u", "data": array(34, 34), "err": null, "mask": null},
    {"label": "reconstructed_psf_u", "filter": "u", "data": array(34, 34), "err": null, "mask": null},
    ...
    {"label": "emline_gflux_6564", "data": array(34, 34), "ivar": array(34, 34), "mask": array(34, 34)},
    {"label": "stellar_vel", "data": array(34, 34), "ivar": array(34, 34), "mask": array(34, 34)},
    {"label": "d4000", "data": array(34, 34), "ivar": array(34, 34), "mask": array(34, 34)}
  ]
}
```

Minimum required fields are:
- provenance: a dictionary of original survey and release of dataset
- observation_id: the unique survey-specific identifier for the observation
- ra: central Right Ascension coordinate
- dec: central Declination coordinate
- spaxel_size: the size of each spaxel in arcsec
- spaxels: a list of individual spaxels in the data cube. Each spaxel has the following minimal fields:

  - flux: the flux array in the spaxel
  - ivar: the uncertainty in the spaxel
  - lsf_sigma: the line-spread-function in the spaxel
  - lambda: the wavelength
  - mask: the pixel quality array
  - x: the x-array element (row) position within the cube
  - y: the y-array element (col) position within the cube

Optional fields can include:

- images: A list of associated images for the data cube.  Can be reconstructed bands or derived analysis maps.  Each image has the following minimal fields:

  - label: a label or name describing the contents of the data image
  - data: the array of data
  - error: optional uncertainty array data
  - mask: optional data quality / mask array data

- extra: an object with survey specific extra data or metadata not strictly necessary but perhaps useful

## Illustrated HuggingFace Dataset generator

The easiest way to add data to the AstroPile is via a [HuggingFace-style dataset generator](https://huggingface.co/docs/datasets/image_dataset#loading-script). Here we'll briefly go over the main parts of the generator, using the [DESI dataloading script](https://github.com/AstroPile/AstroPile_prototype/blob/main/scripts/desi/desi.py) as an example.

First we import the usual suspects (`h5py` and `numpy` for data processing, as well as `itertools` for iterating over series). From HuggingFace we import the `datasets` module, alongside some 'features' that we will later use to define the data type in each column. You may need different columnar features for your dataset, and there is a list [available here](https://huggingface.co/docs/datasets/v2.18.0/en/package_reference/main_classes#datasets.Features).

```python
import datasets
from datasets import Features, Value, Array2D, Sequence
from datasets.data_files import DataFilesPatternsDict
import itertools
import h5py
import numpy as np
```

Optionally in the script preamble we can add some metadata to our dataset, such as a citation pointing to an upstream source, a dataset description, a web link, code licence, and version number. These values will be folded via the `DatasetInfo` method in the `_info` function of our dataloader.

```python
_CITATION = """\
@InProceedings{huggingface:dataset,
title = {A great new dataset},
author={huggingface, Inc.
},
year={2020}
}
"""

_DESCRIPTION = """\
Spectra dataset based on DESI EDR SV3.
"""

_HOMEPAGE = ""

_LICENSE = ""

_VERSION = "0.0.1"
```

We can also add our columnar features in the preamble, to be incorporated into our dataloader later in the script:

```python
_BOOL_FEATURES = [
    "ZWARN"
]

_FLOAT_FEATURES = [
    "Z",
    "ZERR",
    "EBV",
    "FLUX_G",
    "FLUX_R",
    "FLUX_Z",
    "FLUX_IVAR_G",
    "FLUX_IVAR_R",
    "FLUX_IVAR_Z",
    "FIBERFLUX_G",
    "FIBERFLUX_R",
    "FIBERFLUX_Z",
    "FIBERTOTFLUX_G",
    "FIBERTOTFLUX_R",
    "FIBERTOTFLUX_Z",
]
```

Now the fun begins :rocket:. Here we set up a GeneratorBasedBuilder class. We'll go over each part of this class step by step.

```python
class DESI(datasets.GeneratorBasedBuilder):
    """TODO: Short description of my dataset."""

    VERSION = _VERSION
```

The [builder config](https://huggingface.co/docs/datasets/v2.18.0/en/package_reference/builder_classes#datasets.BuilderConfig) defines parameters that are used in the dataset building process, in the AstroPile we are working with `*.hdf5` files so we search for these in our dataset directory with `DataFilesPatternsDict.from_patterns`:

```python
    BUILDER_CONFIGS = [
        datasets.BuilderConfig(
            name="edr_sv3",
            version=VERSION,
            data_files=DataFilesPatternsDict.from_patterns(
                {"train": ["edr_sv3/healpix=*/*.hdf5"]}
            ),
            description="One percent survey of the DESI Early Data Release.",
        ),
    ]

    DEFAULT_CONFIG_NAME = "edr_sv3"

    _spectrum_length = 7781
```

The `_info` function defines the columnar features and other information about our dataset, here we have added comment explanations in-line so that the function flow is obvious.

```python
    @classmethod
    def _info(self):
        # First we add all features common to image datasets.
        # Note that a Sequence requres sub-features so that we can parse it!
        # For the spectrum sequence we have added four float32 Value features
        features = {
            "spectrum": Sequence({
                "flux": Value(dtype="float32"),
                "ivar": Value(dtype="float32"),
                "lsf_sigma":  Value(dtype="float32"),
                "lambda": Value(dtype="float32"),
            }, length=self._spectrum_length)
        }

        # Now we adding all the values from the catalog that we defined earlier
        # in the script, we can add them just like we would do to a normal python
        # dict
        for f in _FLOAT_FEATURES:
            features[f] = Value("float32")

        # Adding all boolean flags
        for f in _BOOL_FEATURES:
            features[f] = Value("bool")

        # Finally we add an object ID for later cross matching and search
        features["object_id"] = Value("string")

        # And we return the above information as a DatasetInfo object,
        # alongside some of the global params we defined in the preamble
        return datasets.DatasetInfo(
            # This is the description that will appear on the datasets page.
            description=_DESCRIPTION,
            # This defines the different columns of the dataset and their types
            features=Features(features),
            # Homepage of the dataset for documentation
            homepage=_HOMEPAGE,
            # License for the dataset if available
            license=_LICENSE,
            # Citation for the dataset
            citation=_CITATION,
        )
```

The [split generator](https://huggingface.co/docs/datasets/image_dataset#download-and-define-the-dataset-splits) function splits our dataset into multiple chunks. Usually this is used for train/test/validation split, but here we define a single 'train' split.

```python
    def _split_generators(self, dl_manager):
        """We handle string, list and dicts in datafiles"""
        if not self.config.data_files:
            raise ValueError(
                f"At least one data file must be specified, but got data_files={self.config.data_files}"
            )
        data_files = dl_manager.download_and_extract(self.config.data_files)
        if isinstance(data_files, (str, list, tuple)):
            files = data_files
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            return [
                datasets.SplitGenerator(
                    name=datasets.Split.TRAIN, gen_kwargs={"files": files}
                )
            ]
        splits = []
        for split_name, files in data_files.items():
            if isinstance(files, str):
                files = [files]
            # Use `dl_manager.iter_files` to skip hidden files in an extracted archive
            files = [dl_manager.iter_files(file) for file in files]
            splits.append(
                datasets.SplitGenerator(name=split_name, gen_kwargs={"files": files})
            )
        return splits
```

Finally we define the example generator generator. This is a generator that yields rows of our dataset according to the structure we defined in our features dict. We want to yield an object ID string with every row of data.

```python
    def _generate_examples(self, files, object_ids=None):
        """Yields examples as (key, example) tuples."""
        for j, file in enumerate(itertools.chain.from_iterable(files)):
            with h5py.File(file, "r") as data:
                if object_ids is not None:
                    keys = object_ids[j]
                else:
                    keys = data["object_id"][:]
                
                # Preparing an index for fast searching through the catalog
                sort_index = np.argsort(data["object_id"][:])
                sorted_ids = data["object_id"][:][sort_index]

                for k in keys:
                    # Extract the indices of requested ids in the catalog 
                    i = sort_index[np.searchsorted(sorted_ids, k)]
                    
                    # Parse spectrum data
                    example = {
                        "spectrum": 
                            {
                                "flux": data['spectrum_flux'][i], 
                                "ivar": data['spectrum_ivar'][i],
                                "lsf_sigma": data['spectrum_lsf_sigma'][i],
                                "lambda": data['spectrum_lambda'][i],
                            }
                    }
                    # Add all other requested features
                    for f in _FLOAT_FEATURES:
                        example[f] = data[f][i].astype("float32")

                    # Add all boolean flags
                    for f in _BOOL_FEATURES:
                        example[f] = not bool(data[f][i])    # if flag is 0, then no problem

                    # Add object_id
                    example["object_id"] = str(data["object_id"][i])

                    yield str(data["object_id"][i]), example
```

To load our newly generated dataset into a downstream script we can again use a HuggingFace tool (`datasets.load_dataset`):

```python
from datasets import load_dataset

print('Preparing DESI dataset') 
dset_desi = load_dataset(
    'AstroPile/desi', # this is the path to the directory containing our dataloading script
     trust_remote_code=True, # we need to enable this so that we can run a custom dataloading script
     num_proc=32, # this is the number of parallel processes
     cache_dir='./hf_cache' # this is the directory where HuggingFace caches the dataset
)
```

Now we are working with a normal HF dataset object, so we can use [all the upstream code as-is](https://huggingface.co/docs/datasets/index).

## Data Architecture

AstroPile will provide utilities to easily generate cross-matched or concatenated versions of these datasets. Below is an example of how the user may generate a cross-matched dataset from the HSC and DECaLS surveys:
```python
import datasets
from astropile import cross_match_datasets

# Load the parent datasets
hsc = datasets.load_building(data_dir='~/data/hsc', config_name='pdr3_dud_22.5')
decals = datasets.load_building(data_dir='~/data/decals', config_name='stein_et_al')

# Build a cross-matched dataset by matching the objects in the two surveys within 1 arcsecond
matched_catalog, matched_dset = cross_match_datasets(hsc, decals, match_radius=1.0)

# The cross-matched dataset can now be saved to disk
matched_dataset.save_to_disk("hsc_meets_decals")

# Push to the HuggingFace Datasets Hub, but be careful not to push large datasets as special
# considerations apply (see https://huggingface.co/docs/hub/repositories-recommendations)
matched_dataset.push_to_hub("my_name/hsc_meets_decals")
```

## Lowel level dataset format

### Limitations of the HuggingFace Paradigm
The native HuggingFace Datasets paradigm has a few limitations for our uses cases where we want to have the flexibility to cross-match and combine datasets easily,
as well as being able to handle datasets that can be in the hundreds of TB in size.
pip install duckdb
  - **Cold storage of the data in Parquet files** Parquet files have nice properties: the are columnars, so that only a subset of columns can be loaded if desired, they are compressed which can reduce the size of data files in some applications, and they can be streamed from online storage which allows fast direct access to some examples of the data. However, the price to pay for compression is that they do not allow for random access, i.e. if I want to retrieve the 1000th image of my file, I actually need to read the 999th first images before I get to the 1000th. This means that if we want to cross-match say 10 spectra with a datasets of postage stamps of size 20 TB, it's going to take many hours because we need to read the entire dataset to find the stamps we are interested in.


  - **Static Parquet then Live Arrrow model** The idea of HF Datasets is that they would be stored online as static Parquet files, but when "prepared" during the `load_dataset` process, they are turned into memory-mapped Arrow tables. At that point it becomes very efficient to search, filter, or in general do random access things in O(1) on the dataset. Arrow is great, it would be very easy to generate sub-datasets and cross-match from the Arrow tables. However, for large-scale datasets of 100s of TB, this means that the entire dataset needs to be copied locally again (essentially, twice the storage space is needed) before I can start subselecting a portion of it.


So, the core issue is that Parquet file do not allow us to do efficient search and filtering, and that we don't necessarily have the time or disk space to turn very large datasets into Arrow format.

### Proposal for AstroPile formatting and storage

The idea is to still use the HuggingFace API as much as possible, but with an alternative to the Parquet storage that retains the ability to process files with random access.

  - **File Format** Instead of Parquet, we propose to rely on traditional HDF5 as the standard for storing large data tables. HDF5 has the nice property of being a stable format, and allow for both column-based access and random access along a column. It makes it very efficient to execute a query on some columns, and then retrieve particular rows of some other columns. We can define an HDF5 DatasetBuilder so that HF Datasets can load HDF5 files just in the same way it would load from Parquet files.


  - **Storage Layout for Efficient Matching and Robust Downloading** Thinking ahead to the case where we have datasets of Billions of entries, and in the hundreds of TB, downloading these datasets will require them to be split in smaller chunks of manageable sizes. If we also want to be able to perform fast searches and cross-matches, it would make sense to have these chunks based on the spatial position of objects on the sky. So, we propose to adopt a simple HealPIX-based splitting of the data using the Hive decomposition strategy. i.e.:
  ```
  root_dataset_folder
    |
    - subset
        |
        - healpix=9999
            |
            - 001-of-002.hdf
            - 002-of-002.hdf
        - healpix=10000
            |
            - 001-of-001.hdf
            ...
  ```

  - **Data Schema for Stored Tables** There is a wide diversity of quantities/features, establishing a common nomenclature for all types of data is not possible. However, we can define a few common columns that are present in all tables and that can be used to perform cross-matches. For instance, we can define a `ra` and `dec` column, and a `healpix` column that can be used to perform spatial searches. For particular types of data, for instance images, spectra, timeseries, we can define a minimal set of columns that are present in all tables of that type. For instance, for images, we can define an `image` column, a `seeing` column, a `pixel_size` column, a `filter` column, etc... These would also have standardized names, to make it easy to access across datasets.
  The tables would be exported in the native official schema for that survey, and the DatasetBuilder would be responsible for applying a set of aliases at reading time to standardize some of the names. The benefit of this approach is that even if we lose the code that exported the dataset, the HDF5 should remain parsable using the official schema.
  Concretely the schema for a dataset exported from HSC would look like:
  ```
  object_id, ra, dec, healpix, forced.g_cmodel_mag, meas.g_cmodel_mag, ..., image, ...
  ```
  i.e. a few names have been standardized, and the rest of the columns are kept as they are in the original survey's official schema.
  The dataset readme file would contain a description of the schema which would be reviewed as part of the dataset review process.
  *Note*: It may also be desirable to define what each survey would consider to be the most general and "go to" say magnitude measurement (for instance forced.g_cmodel_mag - a_g).
