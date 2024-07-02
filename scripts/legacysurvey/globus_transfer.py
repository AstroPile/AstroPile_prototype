import os
import numpy as np
import argparse
import urllib.request
from astropy.table import Table, unique
from globus_sdk import TransferClient, TransferData, NativeAppAuthClient, AccessTokenAuthorizer
from globus_sdk.scopes import TransferScopes

DR10_SOUTH_BRICKS_URL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/survey-bricks-dr10-south.fits.gz"
NERSC_GLOBUS_ENDPOINT = "9d6d994a-6d04-11e5-ba46-22000b92c6ec"


def main(args):
    # Globus endpoint IDs
    source_endpoint_id = NERSC_GLOBUS_ENDPOINT
    destination_endpoint_id = args.destination_endpoint_id
    # User-provided endpoint path
    destination_path = args.destination_path

    # Retrieve the DESI tilepix file
    if not os.path.exists("survey-bricks-dr10-south.fits.gz"):
        print("Downloading bricks info...")
        urllib.request.urlretrieve(DR10_SOUTH_BRICKS_URL, "survey-bricks-dr10-south.fits.gz")
        print("Download complete.")

    # Opening the file and extracting all 
    bricks = Table.read("survey-bricks-dr10-south.fits.gz")
    # Only keep the rows that are primary and with all bands available
    bricks = bricks[bricks["survey_primary"] & (bricks["nexp_g"] > 0) & (bricks["nexp_r"] > 0) & (bricks["nexp_i"] > 0) & (bricks["nexp_z"] > 0)]
    bricks = bricks["brickname"]
            
    # this is the tutorial client ID
    # replace this string with your ID for production use
    CLIENT_ID = "61338d24-54d5-408f-a10d-66c06b59f6d2"
    client = NativeAppAuthClient(CLIENT_ID)

    # Start the login flow
    scope = f"urn:globus:auth:scope:transfer.api.globus.org:all[*https://auth.globus.org/scopes/{NERSC_GLOBUS_ENDPOINT}/data_access]"
    client.oauth2_start_flow(requested_scopes=[TransferScopes.all, scope])
    # client.oauth2_start_flow()

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


    # Splitting the entire data into transfer requests of 10_000 brick each
    # to avoid hitting the maximum number of files per transfer request
    batch_size = 20_000
    for i in range(len(bricks)//batch_size):
        bricks_chunk = bricks[i*batch_size:(i+1)*batch_size]
        if i ==  len(bricks)//batch_size - 1:
            bricks_chunk = bricks[i*batch_size:]

        print("Processing chunk %d of %d bricks..." % (i, len(bricks)//batch_size))

        # Create a transfer data object
        transfer_data = TransferData(transfer_client,
                                    source_endpoint_id,
                                    destination_endpoint_id,
                                    label="Legacy Survey DR10 South Transfer Chunk %d" % i,
                                    sync_level="size")
        
        # Add item for all healpix pixels and all surveys
        n_files = 0
        for brickname in bricks_chunk:
            brick_group = brickname[:3]
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-g.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-g.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-r.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-r.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-i.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-i.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-z.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-image-z.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-g.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-g.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-r.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-r.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-i.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-i.fits.fz")
            n_files += 1
            transfer_data.add_item(f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-z.fits.fz", 
                                destination_path+f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-invvar-z.fits.fz")
            n_files += 1
            transfer_data.add_item(
                f"/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-blobmodel.jpg",
                destination_path
                + f"/dr10/south/coadd/{brick_group}/{brickname}/legacysurvey-{brickname}-maskbits.fits.fz",
            )
            n_files += 1

        print("Submitting transfer request for %d files..." % n_files)
        # Initiate the transfer
        transfer_result = transfer_client.submit_transfer(transfer_data)

        # Get the transfer ID
        transfer_id = transfer_result["task_id"]
        print(f"Transfer submitted with ID: {transfer_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transfer data from DESI to user-provided endpoint.")
    parser.add_argument("destination_endpoint_id", type=str, help="The destination Globus endpoint ID.")
    parser.add_argument("destination_path", type=str, help="The destination path on the endpoint.")
    args = parser.parse_args()
    main(args)
