
import argparse
import pathlib
import urllib

from astropy.io import fits
from globus_sdk import TransferClient, TransferData, NativeAppAuthClient, AccessTokenAuthorizer
from globus_sdk.scopes import TransferScopes


# Main catalog files
MANGA_DRPALL_URL = "https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/drpall-v3_1_1.fits"
MANGA_DAPALL_URL = "https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/dapall-v3_1_1-3.1.0.fits"


# GLOBUS endpoint ID for public SDSS data
SDSS_GLOBUS_ENDPOINT = "f8362eaf-fc40-451c-8c44-50b71ec7f247"



def download_catalog_file(url: str, filepath: str = '.') -> str:
    """ Download a catalog file """
    # Retrieve a catalog file
    name = url.rsplit('/')[-1]
    path = pathlib.Path(filepath) / name

    # return the path if it already exists
    if path.exists():
        print(f'File {path} already exists.')
        return path

    # download the file
    print("Downloading the SDSS MaNGA catalog through HTTPS...")
    downloaded_file, msg = urllib.request.urlretrieve(url, filename=path)
    print("Download complete.")
    return downloaded_file


def _extract_from_drpall(catalog: str, limit: int = None) -> list:
    """ Extract a list of cube file paths from the drpall file """
    base = pathlib.Path('/dr17/manga/spectro/redux/v3_1_1/')

    cubes = []
    with fits.open(catalog) as hdulist:
        n_files = limit or len(hdulist['MANGA'].data)

        for plateifu in hdulist['MANGA'].data['PLATEIFU'][0: n_files]:
            plate, __ = plateifu.split('-')
            path = base / plate / 'stack' / f'manga-{plateifu}-LOGCUBE.fits.gz'
            cubes.append(path)

        return cubes


def _extract_from_dapall(catalog: str, limit: int = None, daptype: str = 'HYB10-MILESHC-MASTARSSP') -> list:
    """ Extract a list of maps file paths from the dapall file """

    base = pathlib.Path('/dr17/manga/spectro/analysis/v3_1_1/3.1.0/')
    daptype = daptype.upper()

    maps = []
    with fits.open(catalog) as hdulist:
        # get the data from the designated DAP-TYPE
        data = hdulist[daptype].data

        # select only files where the DAP was completed
        sub = data[data['DAPDONE'] is True]


        n_files = limit or len(sub)
        for plateifu in sub['PLATEIFU'][0: n_files]:
            plate, ifu = plateifu.split('-')
            path = base / daptype / plate / ifu / f'manga-{plateifu}-MAPS-{daptype}.fits.gz'
            maps.append(path)

        return maps


def build_filelist(product: str, catalog: str, limit: int = None):
    """ Build a list of relevant MaNGA files """
    if product =='cubes':
        files = _extract_from_drpall(catalog, limit)
    elif product == 'maps':
        files = _extract_from_dapall(catalog, limit)
    return files



def transfer_data(destination_endpoint_id: str, destination_filepath: str = '.',
                  client_id: str = None, limit: int = None, product: str = 'cubes'):
    """ Transfer SDSS MaNGA data via Globus """

    # Globus endpoint IDs
    source_endpoint_id = SDSS_GLOBUS_ENDPOINT

    # your globus client id
    client = NativeAppAuthClient(client_id)

    # Start the login flow
    scope = f"urn:globus:auth:scope:transfer.api.globus.org:all[*https://auth.globus.org/scopes/{destination_endpoint_id}/data_access]"
    client.oauth2_start_flow(requested_scopes=[TransferScopes.all, scope])

    # Get the authorization URL
    authorize_url = client.oauth2_get_authorize_url()

    # Print the authorization URL and ask the user to visit it
    print(f'Please go to this URL and login: {authorize_url}')

    # Get the authorization code from the user
    auth_code = input('Please enter the code you get after login here: ').strip()

    # Exchange the authorization code for a token
    token_response = client.oauth2_exchange_code_for_tokens(auth_code)

    globus_transfer_data = token_response.by_resource_server["transfer.api.globus.org"]

    # most specifically, you want these tokens as strings
    TRANSFER_TOKEN = globus_transfer_data["access_token"]

    # Use the access token to create a Globus transfer client
    authorizer = AccessTokenAuthorizer(TRANSFER_TOKEN)
    transfer_client = TransferClient(authorizer=authorizer)

    # Create a transfer data object
    transfer_data = TransferData(transfer_client,
                                source_endpoint_id,
                                destination_endpoint_id,
                                label="SDSS Transfer" ,
                                sync_level="size")

    # setup products
    if product == 'cubes':
        # get catalogs file and build list of manga cubes
        drpall = download_catalog_file(MANGA_DRPALL_URL, filepath=destination_filepath)
        files = build_filelist('cubes', drpall, limit=limit)
    elif product == 'maps':
        # get catalogs file and build list of manga maps
        dapall = download_catalog_file(MANGA_DAPALL_URL, filepath=destination_filepath)
        files = build_filelist('maps', dapall, limit=limit)

    n_files = len(files)

    outpath = pathlib.Path(destination_filepath)
    for file in files:
        outpath = outpath / file
        outpath.parent.mkdir(parents=True, exists_ok=True)

        transfer_data.add_item(file, str(outpath))


    print(f"Submitting transfer request for {n_files} files...")
    # Initiate the transfer
    transfer_result = transfer_client.submit_transfer(transfer_data)

    # Get the transfer ID
    transfer_id = transfer_result["task_id"]
    print(f"Transfer submitted with ID: {transfer_id}")


if __name__ == "__main__":

    # parse cli arguments
    parser = argparse.ArgumentParser(description="Transfer data from public SDSS endpoint to a user-provided endpoint.")
    parser.add_argument("-d", "destination_endpoint_id", type=str, help="The destination Globus endpoint ID.")
    parser.add_argument("-p", "destination_path", type=str, help="The destination path on the endpoint.")
    parser.add_argument("-c", "client_id", type=str, help="Your Globus client id")
    parser.add_argument("-l", "limit", type=int, help="Limit of source files to download.")
    parser.add_argument("-r", "product", type=str, default='cubes', help="The type of MaNGA data product to download")

    args = parser.parse_args()

    # transfer data via globus
    transfer_data(args.destination_endpoint_id, args.destination_path, client_id=args.client_id,
                  limit=args.limit, product=args.product)