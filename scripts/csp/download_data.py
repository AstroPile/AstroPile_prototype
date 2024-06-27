import os
import argparse
import requests
import tarfile

def main(args):

    file_name = "CSP_Photometry_DR3.tgz"
    # Construct the URL to download the file from Zenodo
    url = f"https://csp.obs.carnegiescience.edu/data/{file_name}"

    # Send a GET request to the Zenodo URL to download the file
    response = requests.get(url, stream=True)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Open a file in binary write mode to save the downloaded content
        with open(args.destination_path+file_name, 'wb') as f:
            # Write tar.gz file to destination folder
            f.write(response.content)

        # Unzip tar.gz file
        tar = tarfile.open(args.destination_path+file_name, "r:gz")
        tar.extractall(args.destination_path)
        tar.close()

        # Remove tar.gz file
        os.remove(args.destination_path+file_name)

        # Print a success message if the file is downloaded successfully
        print(f"File downloaded successfully to {args.destination_path}")

        os.rename(f'{args.destination_path}/DR3', f'{args.destination_path}/CSPDR3')
        os.remove(f'{args.destination_path}/CSPDR3/README')
        os.remove(f'{args.destination_path}/CSPDR3/._tab1.dat')
        os.remove(f'{args.destination_path}/CSPDR3/.format_data.py.swp')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transfer data from CSP DR3 to user-provided destination folder.")
    parser.add_argument("destination_path", type=str, help="The destination path to download and unzip the data into.")
    args = parser.parse_args()
    main(args)
