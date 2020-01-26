#!/usr/bin/env python
# -*- coding: utf-8 -*-

### DUST COMPARISON PRE 10^7 ###

import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import csv

### PLOT ###

plt.figure(figsize=(16,5))  
ax=plt.gca()
ax.set_xlim(0,55000)
ax.set_ylim(bottom=0)




### HST FILTERS ###
        
fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/UVUDF_filters.fits'

filter_fits = fits.open(fileName)

#print(catalogue.info())

### UVUDF ###
filters = ['UVUDF-HST_ACS_WFC_F435W',
           'UVUDF-HST_ACS_WFC_F606W',
           'UVUDF-HST_ACS_WFC_F814W',
           'UVUDF-HST_WFC3_IR_F105W',
           'UVUDF-HST_WFC3_IR_F125W',
           'UVUDF-HST_WFC3_IR_F140W',
           'UVUDF-HST_WFC3_IR_F160W']

filter_label = ['F435W',
                'F606W',
                'F814W',
                'F105W',
                'F125W',
                'F140W',
                'F160W']

for i in range(len(filters)):
    
    wl_filter = filter_fits[1].data[filters[i]][0][0]
    data_filter = filter_fits[1].data[filters[i]][0][1]

#    plt.plot(wl_filter, data_filter, label =  filter_label[i])
    plt.fill_between(wl_filter, data_filter, alpha=0.3, label = filter_label[i])

### SPITZER FILTERS ###

spz_wl_1 = []
spz_data_1 = []
spz_wl_2 = []
spz_data_2 = []

with open('/Users/lester/BEAGLE/BEAGLE-general/filters/080924ch1trans_full.txt') as csvfile:
    spam = csv.reader(csvfile, delimiter='\t')
    for row in spam:
        spz_wl_1.append(float(row[0])*10000)
        spz_data_1.append(float(row[1])*2)

with open('/Users/lester/BEAGLE/BEAGLE-general/filters/080924ch2trans_full.txt') as csvfile:
    spam = csv.reader(csvfile, delimiter='\t')
    for row in spam:
        spz_wl_2.append(float(row[0])*10000)
        spz_data_2.append(float(row[1])*2)


plt.fill_between(spz_wl_1, spz_data_1, alpha=0.3, label = 'IRAC 3.6 $\mu m$')
plt.fill_between(spz_wl_2, spz_data_2, alpha=0.3, label = 'IRAC 4.5 $\mu m$')

### SPECTRA ###
    
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Lester_z_shift_2/my_first_test.fits'
data_fits = fits.open(fileName)

redshift = [0, 3, 6, 9]

wl_spec = data_fits[4].data[0][0]
data_spec = data_fits[5].data[0]

sl = 1
norm = max(data_spec)


for i in range(len(redshift)):
    
    ### SPLINE FUN ###
    '''
    x=np.array(wl_spec[::sl])
    y=np.array(data_spec[::sl]/norm)
    x_new = np.linspace(300, 18000, 5000)
    f = interp1d(x, y, kind='quadratic')
    y_smooth=f(x_new)
    plt.plot( (1+redshift[i]) * x_new, y_smooth, label = 'z = %s' % redshift[i])
    '''
    ### ### ###
    
    plt.plot( (1+redshift[i]) * wl_spec[::sl], data_spec[::sl]/norm, label = 'z = %s' % redshift[i])






plt.legend()
plt.show()
























