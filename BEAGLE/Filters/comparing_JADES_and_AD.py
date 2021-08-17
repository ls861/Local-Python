#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 09:19:20 2021

@author: lester
"""

# =============================================================================
# comparing JADES filters with my ASTRODEEP ones
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
matplotlib.rcParams.update({'font.size': 16})


### PLOT ###
fig, ax = plt.subplots(figsize=(16, 6))
#fig.suptitle('ASTRODEEP Filters')


# =============================================================================
# ASTRODEEP
# =============================================================================
wl = []
f = []
fileName = '/Users/lester/BEAGLE/Filter_Files/Ascii_Astrodeep/astrodeep_filters.fits'
data_fits = fits.open(fileName)
filters = ['HST_ACS_WFC_F435W', 'HST_ACS_WFC_F606W', 'HST_ACS_WFC_F814W', 'HST_WFC3_IR_F105W', 'HST_WFC3_IR_F125W', 'HST_WFC3_IR_F140W', 'HST_WFC3_IR_F160W', 'Paranal_HAWKI_Ks', 'Spitzer_IRAC_I1', 'Spitzer_IRAC_I2']
filter_label = ['ACS F435W', 'ACS F606W', 'ACS F814W', 'WFC3 F105W', 'WFC3 F125W', 'WFC3 F140W', 'WFC3 F160W', 'HAWK-I Ks', 'IRAC I1 3.6 $\mu m$', 'IRAC I2 4.5 $\mu m$']

for i in range(len(filters)):
    wl.append(data_fits['TRANSMISSION'].data[filters[i]][0][0])
    f.append(data_fits['TRANSMISSION'].data[filters[i]][0][1])
data_fits.close()

# for i in range(len(filters)):
#     ax.fill_between(wl[i], f[i], alpha=0.7, label=filter_label[i])

for i in range(len(filters)):
    # ax.fill_between(wl[i], f[i], alpha=0.7, label=filter_label[i])
    ax.plot(wl[i], f[i]/max(f[i]), alpha=0.7, label=filter_label[i])
    
# =============================================================================
# JADES
# =============================================================================
wl = []
f = []
fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/JADES_mock_filters_fixed.fits'
data_fits = fits.open(fileName)
filters = ['UVUDF-HST_ACS_WFC_F435W',
           'UVUDF-HST_ACS_WFC_F606W',
           'UVUDF-HST_ACS_WFC_F814W',
           'UVUDF-HST_WFC3_IR_F105W',
           'UVUDF-HST_WFC3_IR_F125W',
           'UVUDF-HST_WFC3_IR_F140W',
           'UVUDF-HST_WFC3_IR_F160W',
           'CANDELS_GS_filter12_IRACch1',
           'CANDELS_GS_filter12_IRACch2']
filter_label = ['F435W',
                'F606W',
                'F814W',
                'F105W',
                'F125W',
                'F140W',
                'F160W',
                'IRAC 3.6 $\mu m$',
                'IRAC 4.5 $\mu m$']

for i in range(len(filters)):
    wl.append(data_fits['TRANSMISSION'].data[filters[i]][0][0])
    f.append(data_fits['TRANSMISSION'].data[filters[i]][0][1])
data_fits.close()

for i in range(len(filters)):
    ax.fill_between(wl[i], f[i], alpha=0.3, label=filter_label[i])
    



# =============================================================================
# REAL Ks filter: https://www.eso.org/sci/facilities/paranal/instruments/hawki/inst/filters/hawki_Knew.dat
# =============================================================================


filename = "/Users/lester/Documents/GitHub/Local-Python/BEAGLE/Filters/hawki_Knew.dat.txt"
ks_data = np.genfromtxt(filename)
wl = ks_data[:,0] * 10
f = ks_data[:,1]
ax.fill_between(wl, f/max(f), alpha=0.3)



# ax.set_xlim(3500,9800) # ACS
# ax.set_xlim(8800,17100) # WFC3
# ax.set_xlim(19000,24000) # Ks
# ax.set_xlim(31000,51000) # IRAC
ax.set_xlim(3500,51000) # ALL


ax.set_ylim(0.0,1.02)
# ax.set_ylim(0.8,1.02)
ax.set_xlabel(r'Wavelength / $\AA$')
ax.set_ylabel('Transmission')
# plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/333_AD_filters.png')
plt.show()

### ###  ###

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)







































