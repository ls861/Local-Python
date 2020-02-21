#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:56:05 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


# =============================================================================
# CALCULATE FLUXES FROM MAGNITUDE FITS FILE
# =============================================================================


filters = ['F090W_NRC_and_OTE_APP', 'F115W_NRC_and_OTE_APP', 'F150W_NRC_and_OTE_APP', 'F200W_NRC_and_OTE_APP', 'F277W_NRC_and_OTE_APP', 'F356W_NRC_and_OTE_APP', 'F444W_NRC_and_OTE_APP']

filter_label = ['F090W', 'F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/fit/mock_catalogue_005_010.fits'
data_fits = fits.open(fileName)

m = []
for i in range(len(filters)):
    m.append(np.array(data_fits['APPARENT MAGNITUDES'].data[filters[i]]))
    
data_fits.close()


f = []
for i in range(len(filters)):
    f.append(10**( (23.9 - m[i]) / 2.5 ))
    
    

# =============================================================================
# FIND ERRORS
# =============================================================================
    
sig_mag5 = np.array([29.4, 29.6, 29.7, 29.8, 29.4, 29.4, 29.1])
    
sig_flux5 = 10**( (23.9 - sig_mag5) / 2.5 )     # microJy
    
sig_flux = sig_flux5 / 5.
    
    

# =============================================================================
# Create ASCII file    
# =============================================================================
    

header_string = '#ID F090W F115W F150W F200W F277W F356W F444W errF090W errF115W errF150W errF200W errF277W errF356W errF444W\n'

file1 = open("JADES_ascii_001.txt","w+")


file1.write(header_string)

for i in range(len(m[0])):
    row = str(i+1)
    
    for j in range(len(filters)):
        row = row + ' ' + str(f[j][i] + np.random.normal(loc=0.0, scale=sig_flux[j]))
    
    for j in range(len(sig_flux)):
        row = row + ' ' + str(sig_flux[j])

    row = row + '\n'    
    file1.write(row)

file1.close()














