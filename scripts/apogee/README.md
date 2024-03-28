# SDSS APOGEE Dataset Collection

This folder contains the scripts and queries used to build the SDSS APOGEE/Milky Way Mapper spectroscopic parent sample, based on 
all infrared spectra at the time of DR17.

## Sample selection

In the current version of the dataset, we retrieve near infrared spectra from SDSS APOGEE. Only the following
cuts are applied:
```
    - TELESCOPE = "apo25m" or "lco25m"
```

## Data preparation

### Downloading data

The script will attempt to download any missing data during the process of compiling the dataset

### Spectra extraction

You can create the parent sample by running the following script:

```bash
python build_parent_sample.py --apogee_data_path=[path to SDSS data] --output_dir=[output directory]
```

e.g. `python build_parent_sample.py /mnt/data_overflow/sdss/apogee/dr17/apogee /mnt/data_overflow/AstroPile/apogee`

### Documentation

- APGEE datamodel https://data.sdss.org/datamodel/
- Bulk data access page: https://www.sdss.org/dr18/data_access/bulk/


