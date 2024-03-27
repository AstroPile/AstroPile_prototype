import logging
import os

from astroquery.utils.tap.core import TapPlus


def query_catalog_via_tap(save_dir:str, max_mag_z:float=20., max_dec:float=-30., dec_interval:float=5., debug:bool=False):
    """
    Query DESI-LS DR10 catalog via TAP and save to parquet shards.
    Many queries to avoid server maxing out (500k rows). Adjust `dec_interval` as needed.

    Could parallelize but it's a bit rude to the server and may not help if server-limited

    Args:
        save_dir (str): Directory to save the catalog shards to.
        max_mag_z (float, optional): Select galaxies brighter than this z magnitude. Defaults to 20.
        min_dec (float, optional): Select galaxies above this declination (degrees). Defaults to -30.
        dec_interval (float, optional): Select galaxies between (min_dec, dec_interval). Degrees. Defaults to 5.
        debug (bool, optional): Run in debug mode (max 1k rows). Defaults to False.

    Returns:
        None
    """

    noirlab = TapPlus(url="https://datalab.noirlab.edu/tap")

    if debug:
        logging.warning('Running in debug mode, max 1k rows, overriding save_dir')
        select_str = 'SELECT TOP 1000'
        save_dir = 'debug'
    else:
        select_str = 'SELECT'

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)


    min_dec = max_dec - dec_interval
    min_dec = max(min_dec, -90)

    logging.info('Downloading with max_mag_z={}, min_dec={}, max_dec={}'.format(max_mag_z, min_dec, max_dec))

    # https://astroquery.readthedocs.io/en/latest/api/astroquery.utils.tap.Tap.html#astroquery.utils.tap.Tap.launch_job
    job = noirlab.launch_job(
        select_str + f"""
        ra, dec, ls_id, brickid, objid, mag_i, mag_r, mag_g, mag_z, dered_mag_i, dered_mag_r, dered_mag_g, dered_mag_z, type, shape_r, fracmasked_r
        FROM ls_dr10.tractor WHERE (mag_z<{max_mag_z:.0f} and "type" !='PSF' and dec>={min_dec:.0f} and dec<{max_dec:.0f})
        """
    )
    r = job.get_results()

    r = r.to_pandas()
    logging.info(f'Sources: {len(r)}')
    assert len(r) < 500000  # server silently refuses to return more rows than this
    save_name = f'ls_dr10_tractor_mag_z_below_{max_mag_z}_dec_{min_dec}_to_{max_dec}.parquet'
    if debug:
        save_name = save_name.replace('.parquet', '_debug.parquet')
    r.to_parquet(os.path.join(save_dir, save_name))


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    save_dir = 'latest'
    for max_dec in range(-90, 90, 5):  # 0 sources at -90 is no problem
        query_catalog_via_tap(
            save_dir=save_dir, max_dec=max_dec, dec_interval=20, debug=True)
    logging.info(f'Finished downloading catalog shards to {os.path.abspath(save_dir)}')
