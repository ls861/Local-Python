#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:31:44 2020

@author: lester
"""

### ### Linear, then delayed exponential SFH ### ###


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits



fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Jan_20_1_035/mock_catalogue.fits'
data_fits = fits.open(fileName)

#info = data_fits.info()
header = data_fits['STAR FORMATION'].header

z   = data_fits['GALAXY PROPERTIES'].data['redshift']
wl  = data_fits['FULL SED WL'].data['wl'][0]                                    # A
f   = data_fits['FULL SED'].data                                                # erg s-1 cm-2 A-1
tau = data_fits['STAR FORMATION BINS'].data['bin_tau']                          # yrs
msa = data_fits['STAR FORMATION'].data['max_stellar_age']                       # yrs
sfr = data_fits['STAR FORMATION'].data['SFR']                                   # yrs

data_fits.close()

tau_size = 14
msa_size = 7

ind = (wl > 500) & (wl < 20000)
xlim = (0, 20000)
ylim = (-21, -17)

time = 1E9*np.linspace(0, 3, len(msa))
time_now = 1.12 * 1E9


plt.plot(time, sfr)
print(sfr)
print(tau)



