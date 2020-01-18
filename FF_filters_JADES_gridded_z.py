#!/usr/bin/env python
# -*- coding: utf-8 -*-

### DUST COMPARISON PRE 10^7 ###

import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import csv

### PLOT ###

plt.figure(figsize=(16,12)) 
plt.xlabel('Wavelength $\lambda$ ($\AA$)', fontsize=14)

ax1=plt.gca()
ax1.set_xlim(3000,55000)
ax1.set_ylim(10**-17, 10**-11)

ax1.set_yscale('log')
ax1.set_ylabel('Energy Density $\lambda F_\lambda$ ($erg\:s^{-1}\:cm^{-2}$)', fontsize=14)


### SPECTRA ###     #5/8 is 10^7, #6/7 is 10^8, #9 is 10^6
    
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Lester_z_shift_grid_neb_8/my_first_test.fits'
data_fits = fits.open(fileName)

wl_spec = data_fits[6].data[0][0]
redshift = data_fits[1].data

sl = 1


#for i in [2]:
for i in range(len(data_fits[7].data)):
        
    data_spec = data_fits[7].data[i]
    # [erg s^-1 cm^-2 A^-1]
    
    wl_spec_z = (1 + redshift[i][0]) * wl_spec[::sl]
    
    plt.plot(wl_spec_z, wl_spec[::sl] * data_spec[::sl], label = 'z = %i' % redshift[i][0], zorder = 2)

plt.legend(loc = 'upper left')


### FILTERS ###
### JADES ### same as UVUDF plus Spitzer

fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/JADES_mock_filters_fixed.fits'
filter_fits = fits.open(fileName)

#print(catalogue.info())

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

filter_centre = [4350, 6060, 8140, 10500, 12500, 14000, 16000, 36000, 45000]    #these are just related to the names of the filters, potentially the mean?

filter_fwhm_centre = [4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 35465.62, 45024.31]

filter_fwhm = [939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 7431.71, 10096.82]

filter_Ew = [823.16, 1771.39, 1886.72, 2371.86, 2674.41, 3569.84, 2750.2, 6836.16, 8649.93]


### Photometric Points ###

data_photo_AB = data_fits[8].data

for i in range(len(data_fits[7].data)):
    
    colour = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
    
    # lambda f_lambda

    data_photo = ( 3631 * (10**(-np.array(data_photo_AB[i]) / 2.5)) ) / (3.34 * (10**4) * np.array(filter_fwhm_centre))
    
    plt.scatter(filter_fwhm_centre, data_photo, marker = 'x', color = 'k', s = 100, zorder = 3)
    plt.errorbar(filter_fwhm_centre, data_photo, xerr=np.array(filter_fwhm)/2, linestyle='None', color = colour[i], zorder = 4)


### PLOT FILTERS ###

ax2 = ax1.twinx()
ax2.set_ylim(0, 1)
ax2.set_ylabel('Transmission', fontsize=14)

for i in range(len(filters)):
    
    wl_filter = filter_fits[1].data[filters[i]][0][0]
    data_filter = filter_fits[1].data[filters[i]][0][1]
    #plt.fill_between(wl_filter, data_filter, alpha=0.3, label = filter_label[i], zorder = 1)
#plt.legend()
    

plt.show()
    

    









