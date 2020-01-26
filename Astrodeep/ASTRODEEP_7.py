#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:33:40 2019

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

def apparent_mag_from_flux(flux):
    '''
    Takes a flux in uJ, and returns apparent magnitude
    '''
    return -2.5 * np.log10((flux * (10**-6))  /  3631)
        

cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)


filter_label = np.array(['B435', 'V606', 'I814', 'Y105', 'J125', 'JH140', 'H160', 'Ks', 'CH1', 'CH2'])

catalog_label =     np.array(['A2744 c',
                              'A2744 p',
                              'M0416 c',
                              'M0416 p',
                              'M0717 c',
                              'M0717 p',
                              'M1149 c',
                              'M1149 p'])


### REALLY BACK TO BASICS INCASE CAT FILE IS CORRUPT ###

M1149_p_z = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_ZPHOT_complete.cat', names=True, dtype=float)
M1149_p_a = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_A.cat', names=True, dtype=float)
M1149_p_b = np.genfromtxt('/Users/lester/BEAGLE/ASTRODEEP/' + 'MACS1149par_B.cat', names=True, dtype=float)

ind = (M1149_p_z['ZBEST'] >= 0) 

print(len(M1149_p_z))
print(len(M1149_p_a))
print(len(M1149_p_b))

M1149_p_z = M1149_p_z[ind]
M1149_p_a = M1149_p_a[ind]
M1149_p_b = M1149_p_b[ind]

print(len(M1149_p_z))
print(len(M1149_p_a))
print(len(M1149_p_b))

##########################
'''
for i in range(len(catalog_label)):
    
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10) = plt.subplots(1, 10, sharey=True, gridspec_kw={'wspace': 0}, figsize=(15, 4))
    fig.suptitle(catalog_label[i])


    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
    
    for j in range(len(filter_label)):
        gt = (cat['b_' + filter_label[j]] > 0.0000001) & (cat[filter_label[j]] > 0) & (cat['field'] == i)
        lt = (cat['b_' + filter_label[j]] > 0.0000001) & (cat[filter_label[j]] < 0) & (cat['field'] == i)
        

        axes[j].scatter(apparent_mag_from_flux(cat['b_' + filter_label[j]][gt]), cat[filter_label[j]][gt], marker='x')
        axes[j].scatter(apparent_mag_from_flux(cat['b_' + filter_label[j]][lt]), abs(cat[filter_label[j]][lt]), marker='x')
        #axes[j].scatter(apparent_mag_from_flux(M1149_p_b[filter_label[j]][M1149_p_b[filter_label[j]] > 0.0000001]), abs(M1149_p_a[filter_label[j]][M1149_p_b[filter_label[j]] > 0.0000001]), marker='x')
        
        #plt.plot([0, 45], [0, 45], marker='.')
        axes[j].set_title(filter_label[j])

plt.show()
##########################   
'''     

# plotting same as before but with magnification taken into account


cat = np.load('cat.npy') # -ve z removed (z=-1), ID > 100000 removed (inc z=99)
filter_uv = np.load('filter_uv.npy')


'''
[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
'''


title = ['Total', 'Clusters', 'Parallels']
tcp = [cat, cat[np.mod(cat['field'], 2) == 0], cat[np.mod(cat['field'], 2) == 1]]


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



# plotting magnification
for i in range(len(tcp)):
    plt.hist(tcp[i]['MAGNIF'], bins=50)
    plt.yscale('log')
    plt.title(title[i])
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











































































