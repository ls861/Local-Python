#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt


cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)

# PARALLELS ONLY

#cat = cat[np.mod(cat['field'], 2) == 1]



#filter_uv = np.load('filter_uv.npy')

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])

catalog_label =     np.array(['A2744 c',
                              'A2744 p',
                              'M0416 c',
                              'M0416 p',
                              'M0717 c',
                              'M0717 p',
                              'M1149 c',
                              'M1149 p'])

'''
[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
'''

'''
for i in range(len(filter_label)):
    cat = cat[cat[filter_label[i]] < 60]
    cat = cat[cat[filter_label[i]] >8]
'''

print(len(cat))


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = cat['ZBEST']
D_l = cosmo.luminosity_distance(z)  # Mpc

'''
ax = np.zeros(len(filter_label))

fig, (ax) = plt.subplots(len(ax), figsize=(8, 20))

#fig.suptitle('PARALLELS')

for i in range(len(filter_label)):

    D_l_sliced = 0    
    catsliced = 0
    
    D_l_sliced = D_l[(cat[filter_label[i]] < 60) & (cat[filter_label[i]] > 8)]
    catsliced = cat[(cat[filter_label[i]] < 60) & (cat[filter_label[i]] > 8)]    
    
    ax[i].hist(    catsliced[filter_label[i]] -  5 * (np.log10(D_l_sliced.value * (10**6)) - 1), bins=30  )
    ax[i].set(ylabel=filter_label[i], yscale='log', xlim=(-40, 0))
'''



for j in range(len(filter_label)):

    D_l_sliced = 0    
    catsliced = 0
    
    D_l_sliced = D_l[(cat[filter_label[j]] < 60) & (cat[filter_label[j]] > 8)]
    catsliced = cat[(cat[filter_label[j]] < 60) & (cat[filter_label[j]] > 8)] 
    
    count = 0
    
    for i in range(len(catsliced)):
        
        if catsliced[filter_label[j]][i] -  5 * (np.log10(D_l_sliced.value[j] * (10**6)) - 1) < -26:
            if D_l_sliced[i].value > 47000:
                print('    ')
                print(filter_label[j])
                print(catsliced[filter_label[j]][i])
                print(D_l_sliced[i].value)
                print(catsliced['field'][i], catsliced['ID'][i])
            count += 1



























