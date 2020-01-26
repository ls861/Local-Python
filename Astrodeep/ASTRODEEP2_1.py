#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

cat = np.load('cat.npy')                    # -ve z removed (z=-1), ID > 100000 removed (inc z=99)
filter_uv = np.load('filter_uv2.npy')       # number of filter, name of filter, apparent mag of filter
uv_m_app_from_cat = filter_uv[2].astype('float64') 

filter_label =      np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
catalog_label =     np.array(['A2744 c', 'A2744 p', 'M0416 c', 'M0416 p', 'M0717 c', 'M0717 p', 'M1149 c', 'M1149 p'])

title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]

'''
[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
'''

# indicies for    5 < app_mag < 80
ind = (uv_m_app_from_cat > 5) & (uv_m_app_from_cat < 80)

# demagnification
uv_m_app_from_cat_demag = uv_m_app_from_cat + (2.5 * np.log10(cat['MAGNIF']))

# Converting to absolute mags using redshift

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = cat['ZBEST']
D_l = (cosmo.luminosity_distance(z)).value * (10**6) #pc
bracket = (1/(1 + z)) * ( (D_l / 10) ** 2 )

uv_m_abs_from_cat       =       uv_m_app_from_cat       - (2.5* np.log10(bracket))
uv_m_abs_from_cat_demag =       uv_m_app_from_cat_demag - (2.5* np.log10(bracket))

### PLOTS USING APP MAGS FROM CATALOG (+ve) ###

# redefined from above to include the sliced "cat"

cat                         = cat[ind]
uv_m_abs_from_cat           = uv_m_abs_from_cat[ind]
uv_m_abs_from_cat_demag     = uv_m_abs_from_cat_demag[ind]

tcp                 = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]
tcp_uv_m_abs        = [uv_m_abs_from_cat, uv_m_abs_from_cat[np.mod(cat['field'], 2) == 0], uv_m_abs_from_cat[np.mod(cat['field'], 2) == 1]]
tcp_uv_m_abs_demag  = [uv_m_abs_from_cat_demag, uv_m_abs_from_cat_demag[np.mod(cat['field'], 2) == 0], uv_m_abs_from_cat_demag[np.mod(cat['field'], 2) == 1]]


count = 0

for k in [tcp_uv_m_abs, tcp_uv_m_abs_demag]:
    
    count += 1
    
    for i in range(len(tcp_uv_m_abs)):
        
        for j in range(10):
            
            plt.hist(k[i][ (tcp[i]['ZBEST'] > j) & (tcp[i]['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2, histtype='step') 
        
        if count == 1:
            plt.title(title[i] + ' - count plotted against Absolute UV Magnitude per redshift bin')
        elif count == 2:
            plt.title(title[i] + ' - DEMAGNIFIED')
        plt.xlabel('Absolute UV Magnitude')
        plt.xlim(-36, 0)
        plt.yscale('log')
        plt.ylabel('Count')
        plt.legend()
        plt.show()
    




















