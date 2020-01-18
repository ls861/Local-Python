#!/usr/bin/env python
# -*- coding: utf-8 -*-

### DUST COMPARISON PRE 10^7 ###

import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import csv

def index(value, list_of_items):
    ''' function that takes a value, finds the closet match in a given list of
    items and then returns the index of where this value is in the list '''
    
    diff = 9999
    index = -1
    
    for i in range(len(list_of_items)):
        if abs(list_of_items[i] - value) < diff:
            index = i
            diff = abs(list_of_items[i] - value)
            
    return index


### PLOT ###

plt.figure(figsize=(16,5))  
ax=plt.gca()
ax.set_xlim(0,55000)
#ax.set_ylim(bottom=0)



### HST FILTERS ###
        
fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/JADES_mock_filters_fixed.fits'
#fileName = '/Users/lester/BEAGLE/BEAGLE-general/filters/UVUDF_filters.fits'
filter_fits = fits.open(fileName)

#print(catalogue.info())

### JADES ### same as UVUDF plus Spitzer
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

filter_centre = [4350,
                6060,
                8140,
                10500,
                12500,
                14000,
                16000,
                36000,
                45000]

for i in range(len(filters)):
    
    wl_filter = filter_fits[1].data[filters[i]][0][0]
    data_filter = filter_fits[1].data[filters[i]][0][1]

#    plt.plot(wl_filter, data_filter, label =  filter_label[i])
    plt.fill_between(wl_filter, data_filter, alpha=0.3, label = filter_label[i])


### SPECTRA ###
    
fileName = '/Users/lester/BEAGLE/BEAGLE-general/results/Lester_z_shift_2/my_first_test.fits'
data_fits = fits.open(fileName)

redshift = [0]
#redshift = [0, 3, 6, 9]

wl_spec = data_fits[4].data[0][0]
data_spec = data_fits[5].data[0]

sl = 1
norm = max(data_spec)


for i in range(len(redshift)):
    
    plt.plot( (1+redshift[i]) * wl_spec[::sl], data_spec[::sl]/norm, label = 'z = %s' % redshift[i])


### Photometric Filters ###

photo = []

for i in range(len(filters)):
    wl_filter = filter_fits[1].data[filters[i]][0][0]
    data_filter = filter_fits[1].data[filters[i]][0][1]
    num = len(wl_filter)
    
    wl_spec = data_fits[4].data[0][0]
    data_spec = data_fits[5].data[0]

    total = 0
    
    for j in range(len(wl_filter)):
        
        if data_filter[j] > 0.2:
            
            ### the if statement was the only way I could get my head around the normalisation issues ###
            total += data_filter[j] * data_spec[index(wl_filter[j], wl_spec)]
    
    avg = total / (num * norm)
    photo.append(avg)
    
print(photo)

plt.scatter(filter_centre, photo)  



plt.legend()
plt.show()
    




















