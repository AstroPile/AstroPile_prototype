import argparse
import os
import tarfile

import requests


def download_tar_file(url, destination_path, file_name):
    response = requests.get(url + file_name, stream=True)
    if response.status_code == 200:
        with open(destination_path + file_name, "wb") as f:
            f.write(response.content)
    else:
        raise ValueError(
            f"Failed to download text file. Status code: {response.status_code}"
        )


def download_text_file(url, destination_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(destination_path, "w") as file:
            file.write(response.text)
    else:
        raise ValueError(
            f"Failed to download text file. Status code: {response.status_code}"
        )


def read_text_file(file_path):
    with open(file_path, "r") as file:
        contents = file.read()
        return contents


def extract_tar_file(destination_path, file_name, *tar_args):
    tar = tarfile.open(os.path.join(destination_path, file_name))
    tar.extractall(destination_path)
    tar.close()
    os.remove(destination_path + file_name)


def cfa_snII_logic(args):
    download_tar_file(
        urls[args.dataset], args.destination_path, file_names[args.dataset]
    )
    extract_tar_file(args.destination_path, file_names[args.dataset])
    os.rename(
        f"{args.destination_path}CFA_SNII_STDSYSTEM_LC.txt",
        f"{args.destination_path}CFA_SNII/STDSYSTEM_LC.txt",
    )
    os.rename(
        f"{args.destination_path}CFA_SNII_NIR_LC.txt",
        f"{args.destination_path}CFA_SNII/NIR_LC.txt",
    )
    os.remove(f"{args.destination_path}CFA_SNII_NATSYSTEM_LC.txt")
    os.remove(f"{args.destination_path}CFA_SNII_STDSYSTEM_STARS.txt")


def cfa_general_logic(args):
    download_text_file(
        urls[args.dataset] + file_names[args.dataset],
        os.path.join(
            args.destination_path, dir_name[args.dataset], file_names[args.dataset]
        ),
    )


file_names = {
    "cfa_snII": "cfa_snII_lightcurvesndstars.june2017.tar",
    "cfa3": "cfa3lightcurves.standardsystem.txt",
    "cfa_SECCSN": "lc.standardsystem.sesn_allphot.dat",
    "cfa4": "cfa4.lc.stdsystem.fi.ascii",
}
urls = {
    "cfa_snII": "https://lweb.cfa.harvard.edu/supernova/fmalcolm2017/",
    "cfa3": "https://lweb.cfa.harvard.edu/supernova/CfA3/",
    "cfa_SECCSN": "https://lweb.cfa.harvard.edu/supernova/",
    "cfa4": "https://lweb.cfa.harvard.edu/supernova/CfA4/",
}
dir_name = {
    "cfa3": "CFA3",
    "cfa_SECCSN": "CFA_SECCSN",
    "cfa4": "CFA4",
    "cfa_snII": "CFA_SNII",
}
survey_specific_logic = {
    "cfa3": cfa_general_logic,
    "cfa_SECCSN": cfa_general_logic,
    "cfa4": cfa_general_logic,
    "cfa_snII": cfa_snII_logic,
}


def main(args):
    os.makedirs(
        os.path.join(args.destination_path, dir_name[args.dataset]), exist_ok=True
    )
    survey_specific_logic[args.dataset](args)
    # Print a success message if the file is downloaded successfully
    print(f"File downloaded successfully to {args.destination_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Transfer data to user-provided destination folder."
    )
    parser.add_argument(
        "destination_path",
        type=str,
        help="The destination path to download and unzip the data into.",
    )
    parser.add_argument(
        "dataset",
        type=str,
        help="The dataset to be downloaded",
    )
    parser.add_argument(
        "-n",
        "--hyphenate-cols",
        nargs="+",
        default=["SPEC_CLASS", "SPEC_CLASS_BROAD", "PARSNIP_PRED", "SUPERRAENN_PRED"],
    )
    args = parser.parse_args()
    try:
        urls[args.dataset]
    except KeyError:
        raise KeyError(f"Dataset argument should be in {urls.keys()}")
    main(args)
