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
# get sfrs from BEAGLE mock and my input to mock
# =============================================================================

# BEAGLE MOCK OUTPUT
fileName = '/Users/lester/Documents/PhD/param_100/mock_100_LE/mock_catalogue_100_LE.fits'
data_fits = fits.open(fileName)
#print(data_fits.info())
sfr = data_fits['STAR FORMATION'].data['SFR']
data_fits.close()

# INPUT TO MOCK
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/LE/mock_MS_parameters_100_LE.fits'
data_fits = fits.open(fileName)
#print(data_fits[1].header)
sfr_MS = data_fits[1].data['sfr']
data_fits.close()

# =============================================================================
# plot 550x550
# =============================================================================

plt.figure(figsize=(10, 10))

plt.xscale('log')
plt.yscale('log')

plt.plot((0, max(sfr_MS)), (0, max(sfr_MS)), zorder=0) # straight line
plt.scatter(sfr_MS, sfr, color='k', marker='x') 
plt.show()

print(np.argmax(abs(sfr_MS-sfr)))
print(sfr_MS[np.argmax(abs(sfr_MS-sfr))])























