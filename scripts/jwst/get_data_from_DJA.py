import numpy as np
import json
import wget
from bs4 import BeautifulSoup
import os
import tarfile
from astropy.table import Table

base_url = 'https://s3.amazonaws.com/grizli-v2/JwstMosaics/v6/index.html'
output_directory = 'cosmos_web/'
try:
    os.chdir(output_directory)
    os.chdir('..')
except:
    print('output directory not found, made the dir.')
    os.mkdir(output_directory)

file = wget.download(base_url)

f = open(file, "r")
file_content = f.read()

soup = BeautifulSoup(file_content, 'html.parser')

cosmos_files = []
for temp in soup.find_all('a'):
    if ('cosweb' in temp['href']) & ('_sci' in temp['href']):
        cosmos_files.append(temp['href'])
    
print('\n')
for i in range(len(cosmos_files)):
    temp = cosmos_files[i].split('_')[0].split('-')
    if temp[-1] != 'clear':
        print('filter: ', temp[-1])
    else:
        print('filter: ', temp[-2])
    print('url: ', cosmos_files[i])

    file = wget.download(cosmos_files[i], out=output_directory)
    
for temp in soup.find_all('a'):
    if ('cosweb' in temp['href']) & ('photoz' in temp['href']):
        print(temp['href'])
        photoz_url = temp['href']

file = wget.download(photoz_url, out=output_directory)

tar = tarfile.open(file)
tar.extractall(path=output_directory)
tar.close()

fnames = os.listdir(output_directory)
for fname in fnames:
    if 'eazypy.zout' in fname:
        phot_table = Table.read(output_directory + fname)
phot_table
