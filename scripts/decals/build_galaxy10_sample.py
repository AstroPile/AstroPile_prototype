import logging
from multiprocessing import Pool
from io import BytesIO
import os
from pathlib import Path
from functools import partial

import numpy as np
import tqdm
import pandas as pd
import requests
import h5py
from astropy.io import fits


from create_images_from_cutout_server import get_download_url




def run_download(df, overwrite, n_processes=4):
    # download_image(train_df.iloc[0])

    pool = Pool(processes=n_processes)
    rows = [row for _, row in df.iterrows()]

    retry_wrapper_partial = partial(retry_wrapper, download_image, overwrite=overwrite)

    # retry_wrapper_partial(rows[0])

    logging.info('Beginning download')
    res = pool.imap_unordered(retry_wrapper_partial, rows)
    list(res)  # iterate
    logging.info('Finished download')


def pack_to_hdf5(df: pd.DataFrame, save_loc: Path):
    # save to hdf5
    with h5py.File(save_loc, 'w') as f:
 
        for col in ['ra', 'dec', 'label']:
            f.create_dataset(col, data=df[col].values)

        f.create_dataset('image_array', (len(df), 3, 224, 224))  # no values yet
        for i, row in tqdm.tqdm(df.iterrows(), total=len(df)):
            with fits.open(row['file_loc']) as hdu:
                f['image_array'][i] = hdu[0].data.astype(np.float32)

    

def retry_wrapper(func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        # logging.error(f'Error in {func.__name__}: {e}')
        logging.error(e)
        return None


def galaxy_to_filename(row):
    return f'{row["h5_index"]}.fits.gz'


def download_image(row, overwrite=False):
    # download these specific images
    save_loc = row['file_loc']
    if not overwrite and os.path.exists(save_loc):
        return

    url = get_download_url(
        row['ra'],
        row['dec'],
        pixel_size=224
    )
    response = requests.get(url, allow_redirects=True)

    # could write directly like so, but prefer .fits.gz instead of .fits
    # with open(f'{base_dir}/images/{row["h5_index"]}.fits', 'wb') as f:
    #     f.write(response.content)

    # read directly into astropy and save as ,fits.gz
    with BytesIO(response.content) as f:
        astropy_fits = fits.open(f)
        astropy_fits.writeto(save_loc, overwrite=True, output_verify='exception')

def main():

    logging.basicConfig(level=logging.INFO)

    base_dir = Path('/Users/user/repos/zoobot-foundation/foundation/datasets/custom_downstream/decals10')
    save_dir = base_dir / 'images'
    debug = True
    n_processes = 4
    overwrite = False
    
    # https://huggingface.co/datasets/mwalmsley/galaxy10_decals/tree/main
    # includes h5_index
    train_df = pd.read_parquet(base_dir / 'decals10_train.parquet')
    test_df = pd.read_parquet(base_dir  / 'decals10_test.parquet')

    # https://astronn.readthedocs.io/en/latest/galaxy10.html

    # clean up some accidental columns
    # these are the original DESI coordinates, TODO check
    train_df = train_df.rename(columns={'ra_x': 'ra', 'dec_x': 'dec'})
    test_df = test_df.rename(columns={'ra_x': 'ra', 'dec_x': 'dec'})
    # print(train_df.head())

    if debug:
        
        # train_df = train_df[:3000]
        # test_df = test_df[:3000]
        train_df = train_df.sample(200, random_state=42).reset_index(drop=True)
        test_df = test_df.sample(100, random_state=42).reset_index(drop=True)

    for (save_loc, df) in [(base_dir/'decals10_train.hdf5', train_df), (base_dir/'decals10_test.hdf5', test_df)]:
        print(f'{save_loc}: {len(df)}')

        df['file_loc'] = save_dir / df.apply(galaxy_to_filename, axis=1)
        run_download(df, overwrite, n_processes)
        pack_to_hdf5(df, save_loc=save_loc)


if __name__ == '__main__':
    main()