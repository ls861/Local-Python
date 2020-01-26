#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)

filter_label =      np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
catalog_label =     np.array(['A2744 c', 'A2744 p', 'M0416 c', 'M0416 p', 'M0717 c', 'M0717 p', 'M1149 c', 'M1149 p'])

'''
[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
'''

filter_uv_file = np.load('filter_uv.npy')
filter_uv = np.array([filter_uv_file, filter_label[filter_uv_file]])
uv_m_app_from_cat = np.zeros(len(filter_uv[1]))

for i in range(len(filter_uv[1])):
    uv_m_app_from_cat[i] = float(cat[filter_uv[1][i]][i])
    
filter_uv = np.array([filter_uv_file, filter_label[filter_uv_file], uv_m_app_from_cat])

print(len(filter_uv[2]))

print(type(uv_m_app_from_cat[0]))

np.save('filter_uv2', filter_uv)




print(type(filter_uv[0]))













'''


# plotting same as before but with magnification taken into account

# I need apparent mags from corrent filter (taking positive only)
# also want just redshifts greater than 1.5 for now.

filter_chosen_a = filter_label[filter_uv] # creates array of filter names 


uv_m_app_from_cat = np.zeros(len(filter_chosen_a))

for i in range(len(filter_chosen_a)):
    uv_m_app_from_cat[i] = cat[filter_chosen_a[i]][i]


plt.hist(uv_m_app_from_cat, bins=500)

# z > 1.5
z15 = cat['ZBEST'] >= 1.5

cat = cat[z15]
uv_m_app_from_cat = uv_m_app_from_cat[z15]
filter_uv = filter_uv[z15]

plt.hist(uv_m_app_from_cat, bins=500)

# app_mag > 5
am_0 = uv_m_app_from_cat > 5

cat = cat[am_0]
uv_m_app_from_cat = uv_m_app_from_cat[am_0]
filter_uv = filter_uv[am_0]

plt.hist(uv_m_app_from_cat, bins=500)

# app_mag > 0
am_80 = uv_m_app_from_cat < 80

cat = cat[am_80]
uv_m_app_from_cat = uv_m_app_from_cat[am_80]
filter_uv = filter_uv[am_80]

plt.hist(uv_m_app_from_cat, bins=500)
plt.show()

plt.hist(uv_m_app_from_cat, bins=50)
plt.show()


print(len(cat))

title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]


# plotting magnification
for i in range(len(tcp)):
    plt.hist(tcp[i]['MAGNIF'], bins=50)
    plt.yscale('log')
    plt.title(title[i])
    plt.show()



uv_m_app_from_cat_demag = uv_m_app_from_cat + (2.5 * np.log10(cat['MAGNIF']))

print(len(uv_m_app_from_cat), len(uv_m_app_from_cat_demag))

### Converting to absolute mags using redshift ###

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = cat['ZBEST']
D_l = cosmo.luminosity_distance(z) 



uv_m_abs_from_cat = uv_m_app_from_cat - ( 5 * np.log10(D_l.value * (10**5)))
uv_m_abs_from_cat_demag = uv_m_app_from_cat_demag - ( 5 * np.log10(D_l.value * (10**5)))



### PLOTS USING APP MAGS FROM CATALOG (+ve) ###

# redefined from above to include the sliced "cat"
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]
tcp_uv_m_abs = [uv_m_abs_from_cat, uv_m_abs_from_cat[np.mod(cat['field'], 2) == 0], uv_m_abs_from_cat[np.mod(cat['field'], 2) == 1]]
tcp_uv_m_abs_demag = [uv_m_abs_from_cat_demag, uv_m_abs_from_cat_demag[np.mod(cat['field'], 2) == 0], uv_m_abs_from_cat_demag[np.mod(cat['field'], 2) == 1]]

for i in range(len(tcp_uv_m_abs)):
    
    for j in range(10):
        
        plt.hist(tcp_uv_m_abs[i][ (tcp[i]['ZBEST'] > j) & (tcp[i]['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2) 
        
    plt.title(title[i] + ' - count plotted against Absolute UV Magnitude per redshift bin')
    plt.xlabel('Absolute UV Magnitude')
    plt.xlim(-36, 0)
    plt.yscale('log')
    plt.ylabel('Count')
    plt.legend()
    plt.show()
    
    for j in range(10):
        
        plt.hist(tcp_uv_m_abs_demag[i][ (tcp[i]['ZBEST'] > j) & (tcp[i]['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2) 
        
    plt.title(title[i] + ' - DEMAGNIFIED')
    plt.xlabel('Absolute UV Magnitude')
    plt.xlim(-36, 0)
    plt.yscale('log')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

### ### ###





plt.hist(cat['ZBEST'][(cat['ZBEST'] > 1.3) & (cat['ZBEST'] < 6)])
plt.show()


print(len(cat['ZBEST'][(cat['ZBEST'] > 1.3) & (cat['ZBEST'] < 6)]))


print(len(cat))
print(len(cat[cat['H160'] < 27.5]))

plt.hist(cat['H160'], bins=50)



print('TEST')

print(len(tcp[0]))
print(len(tcp_uv_m_abs[0]))
print(len(tcp_uv_m_abs_demag[0]))

ind = tcp[0]['ZBEST'] > 5 # ridiculous that this works

cat1 = tcp[0][ind]
filter_uv = filter_uv[ind]
uv_m_abs1 = tcp_uv_m_abs[0][ind]
uv_m_abs_demag1 = tcp_uv_m_abs_demag[0][ind]

cat2 = cat1[uv_m_abs1 < -25]
filter_uv = filter_uv[uv_m_abs1 < -25]

print(len(cat2))



for i in range(len(cat2)):
    print(cat2['field'][i], cat2['ID'][i], cat2['ZBEST'][i], filter_label[filter_uv[i]], cat2[filter_label[filter_uv[i]]][i])



print(len(filter_uv))
print(len(cat2))

'''





















