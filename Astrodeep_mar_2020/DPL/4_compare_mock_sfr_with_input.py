#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:35:19 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# get sfrs from mock with formation z=99, z=999, msa=11 and my original MS
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_catalogue_DPL_001_11.fits'
data_fits = fits.open(fileName)
#print(data_fits.info())
sfr_11 = data_fits['STAR FORMATION'].data['SFR']
data_fits.close()

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_catalogue_DPL_001_13.fits'
data_fits = fits.open(fileName)
#print(data_fits.info())
sfr_13 = data_fits['STAR FORMATION'].data['SFR']
data_fits.close()

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_catalogue_DPL_001_15.fits'
data_fits = fits.open(fileName)
#print(data_fits.info())
sfr_15 = data_fits['STAR FORMATION'].data['SFR']
data_fits.close()

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_004.fits'
data_fits = fits.open(fileName)
#print(data_fits[1].header)
sfr_MS = data_fits[1].data['sfr']
data_fits.close()

# =============================================================================
# plot 550x550
# =============================================================================

plt.figure(figsize=(10, 10))
plt.xlim(0, 550)
plt.ylim(0, 550)

plt.plot((0, 550), (0, 550), zorder=0) # straight line
plt.scatter(sfr_MS, sfr_11, color='k', marker='x', label='formation z=99') 
plt.scatter(sfr_MS, sfr_13, color='r', marker='x', label='formation z=999') 
#plt.scatter(sfr_MS, sfr_15, color='g', marker='x', label='formation z=999') 
plt.legend()
plt.show()

# =============================================================================
# plot 200x200
# =============================================================================

plt.figure(figsize=(10, 10))
plt.xlim(0, 200)
plt.ylim(0, 200)

plt.plot((0, 200), (0, 200), zorder=0) # straight line
plt.scatter(sfr_MS, sfr_11, color='k', marker='x', label='formation z=99') 
plt.scatter(sfr_MS, sfr_13, color='r', marker='x', label='formation z=999') 
#plt.scatter(sfr_MS, sfr_15, color='g', marker='x', label='formation z=999') 
plt.legend()
plt.show()



























