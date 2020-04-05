#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:22:23 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fsize=4
ID = 70

# =============================================================================
# CALCULATE FLUXES FROM MAGNITUDE FITS FILE
# =============================================================================

#column name from apparent mag table in BEAGLE mock (same as left column in config + _APP)
filters = ['HST_ACS_WFC_F435W_APP', 'HST_ACS_WFC_F606W_APP', 'HST_ACS_WFC_F814W_APP', 'HST_WFC3_IR_F105W_APP', 'HST_WFC3_IR_F125W_APP', 'HST_WFC3_IR_F140W_APP', 'HST_WFC3_IR_F160W_APP', 'Paranal_HAWKI_Ks_APP', 'Spitzer_IRAC_I1_APP', 'Spitzer_IRAC_I2_APP']

#column name from ASTRODEEP config file
filter_label = ['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2']
filter_label_plot = ['F435W', 'F606W', 'F814W', 'F105W', 'F125W', 'F140W', 'F160W', 'Ks', 'IRAC 3.6 $\mu m$', 'IRAC 4.5 $\mu m$']
filter_label_old = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_catalogue_DPL_001_15.fits'
data_fits = fits.open(fileName)
#print(data_fits.info())

z_mock = data_fits['GALAXY PROPERTIES'].data['redshift'][ID-1]
wl_spec_mock = data_fits['FULL SED WL'].data['wl'][0]*(1+z_mock)
f_spec_mock = data_fits['FULL SED'].data[ID-1]/(1+z_mock)

appmag_phot_mock = []
for i in range(len(filters)):
    appmag_phot_mock.append(np.array(data_fits['APPARENT MAGNITUDES'].data[filters[i]])[ID-1])
    
data_fits.close()

f_phot_mock = []
for i in range(len(filters)):
    f_phot_mock.append(10**( (23.9 - appmag_phot_mock[i]) / 2.5 ))






# =============================================================================
# PLOT
# =============================================================================

plt.figure(figsize=(2*fsize, fsize))
plt.xlim(0, 50000)
#plt.ylim(1e-21, 1e-15)
#plt.yscale('log')
plt.plot(wl_spec_mock, f_spec_mock*wl_spec_mock)
plt.show()

plt.figure(figsize=(2*fsize, fsize))
plt.xlim(0, 50000)
#plt.ylim(1e-21, 1e-18)
plt.scatter(filter_fwhm_centre, f_phot_mock)
plt.errorbar(filter_fwhm_centre, f_phot_mock, xerr=filter_fwhm/2, linestyle="None")
plt.show()

plt.figure(figsize=(2*fsize, fsize))
plt.xlim(0, 50000)
#plt.ylim(1e-21, 1e-18)
plt.scatter(filter_fwhm_centre, appmag_phot_mock)
plt.errorbar(filter_fwhm_centre, appmag_phot_mock, xerr=filter_fwhm/2, linestyle="None")
plt.show()
















