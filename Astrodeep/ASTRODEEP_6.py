#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(threshold=25, edgeitems=10)
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)
filter_uv = np.load('filter_uv.npy')


'''
[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
'''

filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])
filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])



title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]


plt.figure(figsize=(8, 5))
plt.hist(filter_fwhm_centre[filter_uv] / (1+cat['ZBEST']), bins=30)
plt.show()


plt.figure(figsize=(10, 5))
for i in range(len(filter_label)):
    plt.scatter(cat['ZBEST'][filter_uv == i],   (filter_fwhm_centre[filter_uv] / (1+cat['ZBEST']))[filter_uv == i], label=filter_label[i], s=1, marker='o')
plt.ylim(1300, 4500)
plt.legend()
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(  (filter_fwhm_centre[filter_uv] / (1+cat['ZBEST']))[cat['ZBEST'] > 1.5]  , bins=30)
plt.show()


### ### ### ### ###
### Calculating Absolute mags from chosen filters ###
filter_chosen_a = filter_label[filter_uv] # creates array of filter names 
filter_chosen_b = ['b_' + x for x in filter_chosen_a] # and then appends _b

# again, not obvious how to avoid a loop

uv_flux = np.zeros(len(filter_chosen_a))

for i in range(len(filter_chosen_a)):
    uv_flux[i] = cat[filter_chosen_b[i]][i]
    

cat = cat[uv_flux > 0.0000001]
uv_flux = uv_flux[uv_flux > 0.0000001]

# Histogram of fluxes in chosen filter (>0.0000001) 
logbins = np.geomspace(0.0000001, 25000, 100)

plt.title('Histogram of fluxes in chosen filter ($>$0.0000001)')
plt.xlabel('Flux / $\mu Jy$')
plt.ylabel('Count')
plt.xscale('log')

for i in range(10):
    
    plt.hist(uv_flux[ (cat['ZBEST'] > i) & (cat['ZBEST'] <= i + 1)], bins=logbins, label='%i' % i)
plt.legend()
plt.show()



# calculating absolute mags

uv_m_app = -2.5 * np.log10(  (uv_flux * (10**-6))  /  3631) # calculated apparent mags from flux



# m_app_from_cat is apparent magnitudes from catalog for a comparison with my calculated ones frm the flux
# I'm sure I'm right up to here because all calculated mags match the positive catalog mags

uv_m_app_from_cat = np.zeros(len(cat))

for i in range(len(cat)):
    uv_m_app_from_cat[i] = cat[filter_chosen_a[i]][i]
    
print(uv_m_app)
print(uv_m_app_from_cat)




### Converting fluxes to absolute mags using redshift ###

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = cat['ZBEST']
D_l = cosmo.luminosity_distance(z) 

uv_m_abs = uv_m_app - ( 5 * np.log10(D_l.value * (10**5)))



#removing these high valued mags (assume from 0 flux or something similar, or negative flux?)

cat_from_cat = cat[(uv_m_app_from_cat > 0) & (uv_m_app_from_cat < 80)]
D_l_from_cat = D_l[(uv_m_app_from_cat > 0) & (uv_m_app_from_cat < 80)]
uv_m_app_from_cat = uv_m_app_from_cat[(uv_m_app_from_cat > 0) & (uv_m_app_from_cat < 80)]


uv_m_abs_from_cat = uv_m_app_from_cat[uv_m_app_from_cat > 0] - ( 5 * np.log10(D_l_from_cat.value[uv_m_app_from_cat > 0] * (10**5)))


'''
print(min(uv_m_app_from_cat[uv_m_app_from_cat > 80])) #98.5647
print(max(uv_m_app_from_cat[uv_m_app_from_cat > 80])) #99.0
'''


plt.hist(uv_m_app, bins=np.arange(5, 36, 1))
plt.show()
plt.hist(uv_m_app_from_cat, bins=np.arange(5, 36, 1))
plt.show()


plt.hist(uv_m_abs, bins=np.arange(-35, 6, 1))
plt.show()
plt.hist(uv_m_abs_from_cat, bins=np.arange(-35, 6, 1))
plt.show()

print(len(uv_m_abs))
print(len(uv_m_abs_from_cat))

print(len(cat))
print(len(cat_from_cat))






### PLOTS ###

# redefined from above to include the sliced "cat"
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]
tcp_uv_m_abs = [uv_m_abs, uv_m_abs[np.mod(cat['field'], 2) == 0], uv_m_abs[np.mod(cat['field'], 2) == 1]]


for i in range(len(tcp_uv_m_abs)):
    
    for j in range(10):
        
        plt.hist(tcp_uv_m_abs[i][ (tcp[i]['ZBEST'] > j) & (tcp[i]['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), histtype='step', linewidth=2) 
        
    plt.title(title[i] + ' - count plotted against Absolute UV Magnitude per redshift bin')
    plt.xlabel('Absolute UV Magnitude')
    plt.xlim(-30, 0)
    plt.yscale('log')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

### ### ###


### PLOTS USING APP MAGS FROM CATALOG (+ve) ###

# redefined from above to include the sliced "cat"
tcp = [cat_from_cat, cat_from_cat[np.mod(cat_from_cat['field'], 2) == 0], cat_from_cat[np.mod(cat_from_cat['field'], 2) == 1]]
tcp_uv_m_abs = [uv_m_abs_from_cat, uv_m_abs_from_cat[np.mod(cat_from_cat['field'], 2) == 0], uv_m_abs_from_cat[np.mod(cat_from_cat['field'], 2) == 1]]


for i in range(len(tcp_uv_m_abs)):
    
    for j in range(10):
        
        plt.hist(tcp_uv_m_abs[i][ (tcp[i]['ZBEST'] > j) & (tcp[i]['ZBEST'] <= j+1)], bins=np.arange(-35, 5, 1), label='%i $<$ z $\leq$ %i' % (j, j+1), linewidth=2) 
        
    plt.title(title[i] + ' - count plotted against Absolute UV Magnitude per redshift bin')
    plt.xlabel('Absolute UV Magnitude')
    plt.xlim(-30, 0)
    plt.yscale('log')
    plt.ylabel('Count')
    plt.legend()
    plt.show()

### ### ###














