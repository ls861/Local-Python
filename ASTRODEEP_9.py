#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""



# SANTINI includes 1711 sources in the redshift range 1.3 <= z < 6.


import numpy as np
import matplotlib.pyplot as plt


cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)
filter_uv = np.load('filter_uv.npy')

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


title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]


# CLUSTERS ONLY

filter_uv = filter_uv[np.mod(cat['field'], 2) == 0]
clusters = tcp[1]


print(len(clusters)) # 15296


# For our analysis we applied a magnitude cut of H<27.5.

filter_uv = filter_uv[clusters['H160'] < 27.5]
clusters = clusters[clusters['H160'] < 27.5]


print(len(clusters)) # 9756

# worth seeing what relflag removes:

filter_uv = filter_uv[clusters['RELFLAG'] == 1]
clusters = clusters[clusters['RELFLAG'] == 1]

print(len(clusters)) # 7010

# We have then excluded from our analysis sources at z > 4 whose K and IRAC fluxes have been ignored as their SEDs result highly unconstrained and, as a consequence, their inferred properties are highly unreliable. They amount to 7% of the H-selected sample.

# 1.3 <= z < 6 in the first four HST Frontier Fields

filter_uv = filter_uv[ (clusters['ZBEST'] >= 1.3) & (clusters['ZBEST'] < 11) ]
clusters = clusters[ (clusters['ZBEST'] >= 1.3) & (clusters['ZBEST'] < 11) ]

print(len(clusters)) # 2773
print(len(filter_uv)) # 2773

# After visual inspection of sources with extreme values of the UV slope, we removed sources with beta > 1 or beta >= −3.5, mostly caused by noisy photometry (~8% of the sample in the redshift range analysed).

# We adopted the median magnification

# removed bunch of items due to mass limit 

# After excluding 17 objects at z > 3 which have been a-posteriori visually inspected

# SANTINI includes 1711 sources in the redshift range 1.3 <= z < 6.




### PLOTS ###


filter_chosen_a = filter_label[filter_uv] # creates array of filter names 


uv_m_app_from_cat = np.zeros(len(filter_chosen_a))

for i in range(len(filter_chosen_a)):
    uv_m_app_from_cat[i] = cat[filter_chosen_a[i]][i]



title = 'Clusters'
uv_m_app_from_cat_demag = uv_m_app_from_cat + (2.5 * np.log10(clusters['MAGNIF']))

### Converting to absolute mags using redshift ###


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = clusters['ZBEST']
D_l = cosmo.luminosity_distance(z)  # Mpc


uv_m_abs_from_cat = uv_m_app_from_cat -  5 * (np.log10(D_l.value * (10**6)) - 1) 
uv_m_abs_from_cat_demag = uv_m_app_from_cat_demag - 5 * (np.log10(D_l.value * (10**6)) - 1) 


print(uv_m_app_from_cat)
print(cosmo.luminosity_distance(3).value )

plt.hist(uv_m_app_from_cat, bins = 50)
plt.show()

plt.hist(uv_m_app_from_cat_demag, bins = 50)
plt.show()

plt.hist(uv_m_abs_from_cat, bins = 50)
plt.show()

plt.hist(uv_m_abs_from_cat_demag, bins = 50)
plt.show()



### PLOTS USING APP MAGS FROM CATALOG (+ve) ###

# redefined from above to include the sliced "cat"

tcp_uv_m_abs = uv_m_abs_from_cat
tcp_uv_m_abs_demag = uv_m_abs_from_cat_demag


    
for j in range(10):
    
    plt.hist(tcp_uv_m_abs[ (clusters['ZBEST'] > j) & (clusters['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2) 
    
plt.title(title + ' - count plotted against Absolute UV Magnitude per redshift bin')
plt.xlabel('Absolute UV Magnitude')
plt.xlim(-36, 0)
plt.yscale('log')
plt.ylabel('Count')
plt.legend()
plt.show()

for j in range(10):
    
    plt.hist(tcp_uv_m_abs_demag[ (clusters['ZBEST'] > j) & (clusters['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2) 
    
    
plt.title(title + ' - DEMAGNIFIED')
plt.xlabel('Absolute UV Magnitude')
plt.xlim(-36, 0)
plt.yscale('log')
plt.ylabel('Count')
plt.legend()
plt.show()

### ### ###



print(len(tcp_uv_m_abs_demag))
print(len(clusters))

test = clusters[(tcp_uv_m_abs_demag < -25) & (clusters['ZBEST'] > 5)]
print(test['field'], test['ID'])

print(len(test))





































