
import argparse
import os

from sdss_access import Access


def transfer_data(limit: bool = None, daptype: str = 'HYB10-MILESHC-MASTARSSP',
                  test: bool = False) -> None:
    """ Download SDSS MaNGA files

    Use parallelized rsync streams to download SDSS data
    with sdss_access.  This will download all files unless the
    limit boolean flag is set.  The limit boolean flag will only
    download a subset of data for testing, for a single SDSS plate.

    For the MaNGA DAP analysis types, see
    "https://www.sdss4.org/dr17/manga/manga-data/data-model/#DataAnalysisPipelineOutput"

    Parameters
    ----------
    limit : bool, optional
        Flag to limit the data downloaded, by default None
    daptype : str, optional
        The SDSS MaNGA DAP analysis type, by default "HYB10-MILESHC-MASTARSSP"
    test : bool, optional
        Flag to test without actually downloading
    """
    access = Access(release='DR17')
    access.remote()

    # get catalog files
    access.add('drpall', drpver='v3_1_1')
    access.add('dapall', drpver='v3_1_1', dapver='3.1.0')

    # get all cube / maps files
    if limit:
        access.add('mangacube', drpver='v3_1_1', wave='LOG', plate='8485', ifu='*')
        access.add('mangadap', plate='8485', ifu='*', drpver='v3_1_1', dapver='3.1.0', mode='MAPS', daptype=daptype)
    else:
        access.add('mangacube', drpver='v3_1_1', wave='LOG', plate='*', ifu='*')
        access.add('mangadap', plate='*', ifu='*', drpver='v3_1_1', dapver='3.1.0', mode='MAPS', daptype=daptype)

    access.set_stream()

    print(f'Downloading {len(access.stream.task)} files.')
    if not test:
        access.commit()


if __name__ == "__main__":

    # parse cli arguments
    parser = argparse.ArgumentParser(description="Transfer SDSS data with parallelized rsync using sdss_access")
    parser.add_argument("-p", "destination_path", type=str, default='.', help="The destination directory")
    parser.add_argument("-l", "limit", action='store_true', default=False, help="Flag to download a subset of files.")

    args = parser.parse_args()

    # set the top level SDSS download directory to the destination path
    # this envvar is needed for sdss_access.  All file downloads paths are relative to this
    # directory and retain the directory structure of the SDSS Science Archive Server (SAS)
    os.environ["SAS_BASE_DIR"] = args.destination_path

    # transfer data via globus
    transfer_data(limit=args.limit)
