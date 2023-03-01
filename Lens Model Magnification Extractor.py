#!/usr/bin/env python
# coding: utf-8

# # Lens Model Magnification Extractor

# This script returns lens model magnification estimates from three Hubble Treasury programs: CLASH, Frontier Fields, and RELICS. Inputs are RA, Dec, z of the strongly lensed galaxy. The script determines the nearest galaxy cluster and finds all lens models on MAST, reading the files in place on the server (without downloading to your computer), and then outputs the magnification estimates. In addition to the "best" estimate from each model, the user may optionally (not by default) choose to also extract the full range of ~100 estimates from each lens model to probe the uncertainties.

# WARNING: This will download model FITS files to your Astropy cache (a hidden directory with cryptic filenames). If you analyze many clusters, this could consume hundreds of GB of disk space. A better implementation of this code would download the models to a non-hidden directory preserving filenames and directory structures...

# The code is the same as used by the Frontier Fields lens model interactive web tool:    
# https://archive.stsci.edu/prepds/frontier/lensmodels/webtool

# The penultimate command lists the available lens models. Users should take care to select desired models, including more recent versions with better constraints (e.g., v2 instead of v1). Accessing each online file takes some time. The Frontier Fields have many lens models contributed by various groups; it would take a long time to access them all. The FF webtool runs more quickly because the files are on the same server.

# Questions / feedback / suggestions: Dan Coe dcoe@stsci.edu

# ### CLASH
# https://www.stsci.edu/~postman/CLASH/    
# https://archive.stsci.edu/prepds/clash/
# 
# ### Frontier Fields
# https://outerspace.stsci.edu/display/HPR/HST+Frontier+Fields    
# https://archive.stsci.edu/prepds/frontier/
# 
# ### RELICS
# https://relics.stsci.edu    
# https://archive.stsci.edu/prepds/relics/

from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

import os
import numpy as np

def roundint(x):
    return np.round(x).astype(int)


# Read online directories and return files with names including search terms

import requests
from bs4 import BeautifulSoup

def process_string_list(string_list):
    if type(string_list) == str:
        if ' ' in string_list:
            string_list = string_list.split()
        else:
            string_list = [string_list]
    return string_list

def online_files(url, search_strings=[], exclude_strings=[]):
    '''Return directory list of files at url
    pruned such that filenames contain search_strings and not exclude_strings'''
    soup = BeautifulSoup(requests.get(url).text)
    image_files = [tag.text for tag in soup.find_all('a')]
    search_strings = process_string_list(search_strings)
    for search_string in search_strings:
        image_files = [file for file in image_files if search_string in file]
    exclude_strings = process_string_list(exclude_strings)
    for exclude_string in exclude_strings:
        image_files = [file for file in image_files if exclude_string not in file]
    image_files = [os.path.join(url, file) for file in image_files]
    return image_files

def online_file(url, search_strings=[], exclude_strings=[]):
    '''Return single file from url directory with filename including search_strings and not exclude_strings'''
    image_files = online_files(url, search_strings, exclude_strings)
    if len(image_files) == 1:
        return image_files[0]
    elif len(image_files) == 1:
        raise ValueError('No matches')
    else:
        raise ValueError('More than one match')


import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def Dds_Ds(zl, zs):
    Dds = cosmo.angular_diameter_distance_z1z2(zl, zs)
    Ds  = cosmo.angular_diameter_distance_z1z2(0 , zs)
    return (Dds / Ds).value


def extract_magnification(kappa_file, gamma_file, RA, Dec, z_lens, z_source, verbose=True):
    if verbose:
        print('Loading %s ...' % kappa_file)
        
    hdu = pyfits.open(kappa_file, memmap=1)
    kappa = data = hdu[0].data

    header = hdu[0].header
    wcs = pywcs.WCS(header)
    sky = np.array([[RA, Dec]])
    x, y = wcs.wcs_world2pix(sky, 1)[0]  # with (1,1) at the origin
    
    x = roundint(x)
    y = roundint(y)
    
    ny, nx = data.shape
    if (x < 0) or (x > nx-1) or (y < 0) or (y > ny-1):
        return None
    
    if verbose:
        print('Loading %s ...' % gamma_file)
        
    gamma = pyfits.open(gamma_file, memmap=1)[0].data

    kappa = kappa[y-1, x-1]
    gamma = gamma[y-1, x-1]

    Dds_Ds1 = Dds_Ds(z_lens, z_source)
    kappa = Dds_Ds1 * kappa
    gamma = Dds_Ds1 * gamma
    
    magnif = 1 / ((1 - kappa)**2 - gamma**2)
    magnif = abs(magnif)
    
    return magnif


from astropy.io import ascii
clusters = ascii.read('http://www.stsci.edu/~dcoe/clusters.txt')


def find_closest_cluster(target_coords, distance_match=5*u.arcmin):
    for i in range(len(clusters)):
        cluster_coords = SkyCoord(clusters['RA'][i], clusters['Dec'][i], unit=(u.hourangle, u.deg))
        distance = target_coords.separation(cluster_coords)
        if distance < distance_match:
            cluster = clusters['name'][i]
            z_lens = float(clusters['z'][i])
            program = clusters['program'][i]
            return cluster, z_lens, program
    return

def extract_magnifications(model_directory, output_file=None, confidence_percentile=None):
    print('Directory', model_directory)
    model_name = model_directory.replace(model_top_directory, '')
    if model_name[-1] == '/':
        model_name = model_name[:-1]
        
    model_name = model_name.replace('/', '_')
    
    kappa_file = online_file(model_directory, '_kappa.fits')
    gamma_file = online_file(model_directory, '_gamma.fits')
    magnif_best  = extract_magnification(kappa_file, gamma_file, ra_deg, dec_deg, z_lens, z_source)
    print(ra_deg, dec_deg, z_source, cluster, z_lens, model_name, magnif_best)

    outline = ''
    outline += '  %s' % id.ljust(20)
    outline += '  %11.7f' % ra_deg
    outline += '  % 11.7f' % dec_deg
    outline += '  %6.4f' % z_source
    outline += '  %s' % cluster.ljust(15)
    outline += '  %6.4f' % z_lens

    outline += '  %s' % model_name.ljust(20)
    outline += ' %5.2f' % magnif_best

    magnifications = []

    if confidence_percentile:
        # Extract estimates from range of ~100 models
        range_directory = os.path.join(model_directory, 'range')
        print(range_directory)
        imodel = 0
        while 1:
            try:
                kappa_file = online_file(range_directory, 'map%03d_kappa.fits' % imodel)
                gamma_file = online_file(range_directory, 'map%03d_gamma.fits' % imodel)
                magnif  = extract_magnification(kappa_file, gamma_file, ra_deg, dec_deg, z_lens, z_source)
                print('map%03d' % imodel, magnif)
                magnifications.append(magnif)
                imodel += 1
            except:  # until done
                print('NOT FOUND')
                break
    
    if len(magnifications):
        # confidence_percentile = 68.3  # percent of models nearest median
        sigma1 = (100 - confidence_percentile) / 2
        percentiles = lo, med, hi = np.percentile(magnifications, (sigma1, 50, 100-sigma1))
        outline += ' %5.2f' % magnif_med
        outline += ' %5.2f' % magnif_lo
        outline += ' %5.2f' % magnif_hi

        for magnif in np.sort(magnifications):
            outline += ' %5.2f' % magnif

    outline += '\n'

    if output_file:
        fout = open(output_file, 'a')
        fout.write(outline)
        fout.close()
    
    return outline

# Example inputs

# Frontier Fields (takes priority as most recent) / CLASH 
id = 'MACS1149-JD'
RA  =  '11 49 33.584'
Dec = '22 24 45.78'
z_source = 9.11
target_coords = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
ra_deg = target_coords.ra.value
dec_deg = target_coords.dec.value


# RELICS
id = 'SPT0615-JD'
RA  =  '06 15 55.03'
Dec = '-57 46 19.56'
z_source = 10
target_coords = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
ra_deg = target_coords.ra.value
dec_deg = target_coords.dec.value


# CLASH
id = 'MACS0647-JD1'
RA  =  '06:47:55.731'.replace(':', ' ')
Dec = '+70:14:35.76'.replace(':', ' ')
z_source = 11
target_coords = SkyCoord(RA, Dec, unit=(u.hourangle, u.deg))
ra_deg = target_coords.ra.value
dec_deg = target_coords.dec.value
ra_deg, dec_deg


cluster, z_lens, program = find_closest_cluster(target_coords)
print(cluster, z_lens, program)


def extract_model_directories(directory):
    model_files = online_files(directory, '.fits')
    if len(model_files):
        model_directories.append(directory)

    # Models in subdirectories (e.g., v2)
    for subdirectory in online_files(directory, '/', 'range'):
        extract_model_directories(subdirectory)

model_top_directory = 'https://archive.stsci.edu/missions/hlsp/%s/%s/models/' % (program, cluster)
model_directories = []
extract_model_directories(model_top_directory)
model_directories

headline = 'id                    ra            dec         z_source old_mu cluster         z_lens  lens_model            best  med   lo    hi '
for imodel in range(100):
    headline += '  m%03d' % imodel

output_file = 'test_magnifications.cat'

fout = open(output_file, 'w')
fout.write(headline+'\n')
fout.close()

for model_directory in model_directories:
    extract_magnifications(model_directory, output_file)


