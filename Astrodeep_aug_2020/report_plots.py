#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 12:05:16 2020

@author: lester
"""


import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import pickle
from numpy import errstate,isneginf
import sys



AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
print(AD.dtype.names)

idx_c = np.isin(AD['field'], [0.0, 2.0, 4.0, 6.0])
idx_p = np.isin(AD['field'], [1.0, 3.0, 5.0, 7.0])

AD_c = AD[idx_c]
AD_p = AD[idx_p]

print(AD_p['field'][0])


plt.hist(AD['MAGNIF'], bins=50, range=(1, 80), log=True)
plt.show()
plt.hist(AD_c['MAGNIF'], bins=50, range=(1, 80), log=True)
plt.show()
plt.hist(AD_p['MAGNIF'], bins=50, range=(1, 2), log=True)
plt.show()



# =============================================================================
# -ve mags
# =============================================================================

print(sum(AD['H160']<0))
print(sum(AD['b_H160']<0))

print(sum(AD['H160']<-10))
print(sum(AD['b_H160']<-10))

print(sum(AD['H160']==0))
print(sum(AD['b_H160']==0))

print(sum(AD['H160']==99))
print(sum(AD['b_H160']==99))


plt.hist(AD['H160'], bins=60, range=(-40,-10))
plt.show()

plt.hist(AD['H160'], bins=60, range=(-1,1))
plt.show()

plt.hist(AD['H160'], bins=60, range=(10,40))
plt.show()

plt.hist(AD['H160'], bins=60, range=(97,100))
plt.show()

plt.hist(AD['b_H160'], bins=60, range=(-5,30))
plt.show()

plt.hist(AD['b_H160'], bins=60, range=(-1,1))
plt.show()

plt.hist(AD['b_H160'], bins=60, range=(10,40))
plt.show()

plt.hist(AD['b_H160'], bins=60, range=(97,100))
plt.show()


# =============================================================================
# redshift BEAGLE vs AD
# =============================================================================


fields = ['A2744_c', 'A2744_p', 'M0416_c', 'M0416_p', 'M0717_c', 'M0717_p', 'M1149_c', 'M1149_p']

for i, field in enumerate(fields):
    
    #####
    
    summCat = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(i)

    data_fits = fits.open(summCat)
    ids = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    z_BEAGLE = data_fits['GALAXY PROPERTIES'].data['redshift_median']
#    data = {}
#    data['mTot'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_tot_median'], dtype=np.float64))}
#    data['mStar'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_star_median'], dtype=np.float64))}
#    temp_sfr = np.asarray(data_fits['STAR FORMATION'].data['SFR_median'], dtype=np.float64)
#    data['sfr'] = {'value':np.log10( np.where(temp_sfr==0, 1e-30, temp_sfr) )}
    data_fits.close()
    
    if i == 0:
        z_BEAGLE_arr = z_BEAGLE 
    else:
        z_BEAGLE_arr = np.concatenate((z_BEAGLE_arr, z_BEAGLE))
        
        
    #####
        
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/data/astrodeep_{}_subset_RF1_001.fits'.format(field)
    data_fits = fits.open(fileName)
#    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
#    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
#    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    z_AD = np.array(data_fits[1].data['ZBEST'])
    data_fits.close()
    
    print('LESTER')
    print(len(z_BEAGLE), len(z_AD), len(ids))
    print(len(z_AD[ids-1]))
    
    if field == 'A2744_c':
#        id_input_arr = id_input
#        field_original_arr = field_original
#        id_original_arr = id_original
        z_AD_arr = z_AD[ids-1]
        
    else:
#        id_input_arr = np.concatenate((id_input_arr, id_input))
#        field_original_arr = np.concatenate((field_original_arr, field_original))
#        id_original_arr = np.concatenate((id_original_arr, id_original))
        z_AD_arr = np.concatenate((z_AD_arr, z_AD[ids-1]))

    
print(len(z_AD_arr))
print(len(z_BEAGLE_arr))

plt.scatter(z_BEAGLE_arr, z_AD_arr)
plt.show()

plt.hist2d(z_BEAGLE_arr, z_AD_arr, bins=20)
plt.colorbar()
plt.show()



