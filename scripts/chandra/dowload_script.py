# Astropile - Rafael Martinez-Galarza
# This code dowloads data files from the Chandra Source Catalog
# using the CLI protocol:
# https://cxc.cfa.harvard.edu/ciao/threads/csccli/
# For each X-ray detection (each one associated to a source name 
# and RA and DEC, it brings down the folowing files: 
# 1) Event file - this is the file containing all photon recordings
# 2) Light curve - binned light curves from Gregory-Loredo alg.
# 3) Spectral files - PHA (signal), and ARF and RMF (response)
# 4) Images - fits files containing the counts binned spatially
# 5) PSF - the computed PSF at the position of the detection
# The data downloaded (IDs compiled in auxiliary file obsid_obi_regid.txt)
# correspond to X-ray sources with enough counts to allow for spectral
# fitting (and a minimum of 50 counts), that have S/N of at least 5, 
# and that are at most 5 arcmin off-axis in the field of view.
# Also, all events in the broadband (0.5keV-7.0keV) are included

import requests

# We use the current ('cur') of the Chandra Source Catalog Database
def retrieve(url, packageset, idx):
    response = requests.get(url, params={
        'version': 'cur',
        'packageset': packageset
    })
    with open(f'package.{idx}.tar', 'wb') as output:
        output.write(response.content)

# This is the url for retrieval to data at the CfA
url = 'http://cda.cfa.harvard.edu/csccli/retrieve'

# We will download the data in packages of 50 detections each
number_of_identifiers_per_request = 50

packageset = ''
number_of_identifiers = 0

# The file below contains the list of detection IDs
separator = ''
with open('obsid_obi_regid.txt', 'r') as input:
    # read header line and ignore it
    input.readline()
    
    while True:
        line = input.readline()
        if '' == line:
            break
        line = line.rstrip()

        number_of_identifiers += 1
        
        # Here we specify which datatypes to download
        packageset += separator + line + '/lightcurve/b'
        separator = ','
        packageset += separator + line + '/spectrum/b'
        separator = ','
        packageset += separator + line + '/rmf/b'
        separator = ','
        packageset += separator + line + '/arf/b'
        separator = ','
        packageset += separator + line + '/regimg/b'
        separator = ','
        packageset += separator + line + '/regevt3/b'
        separator = ','
        packageset += separator + line + '/psf/b'
        
        if 0 == number_of_identifiers % number_of_identifiers_per_request:
            retrieve(url, packageset, int(number_of_identifiers / number_of_identifiers_per_request))

            packageset = ''
            separator = ''
            
            # Print progress (about 55,000 detections in total) with
            # the current set thresholds of S/N, etc.
            print(number_of_identifiers)
        

    if 0 != number_of_identifiers % number_of_identifiers_per_request:
        retrieve(url, packageset, int(number_of_identifiers / number_of_identifiers_per_request)+1)
