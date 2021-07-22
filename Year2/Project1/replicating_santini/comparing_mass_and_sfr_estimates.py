#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 16:56:46 2020

@author: lester
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle


AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/astrodeep/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)

sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'

mag_GMM = np.array([np.log10(AD['MAGNIF'])]*3).transpose() # this is genius

with np.errstate(divide='ignore', invalid='ignore'):
    data = {    'field_AD':             AD['field'],
                'id_AD':                AD['ID'], 
                'mag_AD':               np.log10(AD['MAGNIF']), 
                'redshift_AD':          AD['ZBEST'], 
                'mass_AD':              np.log10(AD['MSTAR']*1e9) - np.log10(AD['MAGNIF']), 
                'mass_AD_neb':          np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']), 
                'sfr_AD':               np.log10(AD['SFR']) - np.log10(AD['MAGNIF']), 
                'sfr_AD_neb':           np.log10(AD['SFR_NEB']) - np.log10(AD['MAGNIF']), 
                'relflag_AD':           AD['RELFLAG'],
                'id_BEAGLE':            np.load(sbf+'id_BEAGLE.npy'), 
                'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_tot':      np.log10(np.load(sbf+'mass_BEAGLE_tot.npy')) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_stellar':  np.log10(np.load(sbf+'mass_BEAGLE_stellar.npy')) - np.log10(AD['MAGNIF']), 
                'sfr_BEAGLE_instant':   np.load(sbf+'sfr_BEAGLE_instant.npy') - np.log10(AD['MAGNIF']), 
                'redshift_BEAGLE':      np.load(sbf+'redshift_BEAGLE.npy'), 
                'tau_BEAGLE':           np.load(sbf+'tau_BEAGLE.npy'), 
                'tauv_BEAGLE':          np.load(sbf+'tauv_BEAGLE.npy'), 
                'msa_BEAGLE':           np.load(sbf+'msa_BEAGLE.npy'), 
                'metallicity_BEAGLE':   np.load(sbf+'metallicity_BEAGLE.npy'),
                'min_chi2_BEAGLE':      np.load(sbf+'min_chi2_BEAGLE.npy'),
                
                'id_GMM':               np.load(sbf+'id_GMM.npy'),
                'x_GMM':                np.load(sbf+'x_GMM.npy') - mag_GMM,
                'y_GMM':                np.load(sbf+'y_GMM.npy') - mag_GMM,
                'xsig_GMM':             np.load(sbf+'xsig_GMM.npy'),
                'ysig_GMM':             np.load(sbf+'ysig_GMM.npy'),
                'xycov_GMM':            np.load(sbf+'xycov_GMM.npy'),
                'amp_GMM':              np.load(sbf+'amp_GMM.npy')

                }
    

idx1 = (AD['field']%1.0==0.0) # all
idx2 = (AD['H160']<27.5)
idx = np.logical_and(idx1,idx2)

zLow = 1.5
zHigh = 1.8
wLim = (zHigh - zLow) / 2.0
zLim = zLow + wLim
idx3 = (abs(AD['ZBEST']-zLim) < wLim)
idx = np.logical_and(idx,idx3)

idx5 = (AD['RELFLAG']==1.0)
idx = np.logical_and(idx,idx5)

idx7 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
idx = np.logical_and(idx,idx7)

for key in data.keys():
    data[key] = data[key][idx]

    
#print(len(data['mass_AD']))

#%%

# =============================================================================
# Santini
# =============================================================================
    
san_bin = 0
    
# logSFR = alpha log(M / M_9.7) + beta
z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

alpha_san = B_san - 9.7*A_san
beta_san = A_san
# CHECK THIS
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_err_san = A_err_san

print(alpha_san)
print(beta_san)
# =============================================================================
# THINGS TO COMPARE
# =============================================================================


#start_mass = data['mass_AD']
start_mass = data['mass_AD_neb']
#start_mass = data['mass_BEAGLE_tot']
#start_mass = data['mass_BEAGLE_stellar']

#start_sfr = data['sfr_AD']
#start_sfr = data['sfr_AD_neb']
#start_sfr = data['sfr_BEAGLE_instant']
start_sfr = data['sfr_SAN']

#end_mass = data['mass_AD']
#end_mass = data['mass_AD_neb']
#end_mass = data['mass_BEAGLE_tot']
end_mass = data['mass_BEAGLE_stellar']

#end_sfr = data['sfr_AD']
#end_sfr = data['sfr_AD_neb']
#end_sfr = data['sfr_BEAGLE_instant']
end_sfr = data['sfr_SAN']


#%%

# =============================================================================
# some cuts
# =============================================================================

min_BEAGLE_sfr_cut = -4.0
max_BEAGLE_sfr_cut = 4.0

min_BEAGLE_mass_cut = 5.0
max_BEAGLE_mass_cut = 12.0

idx_cut = [(start_mass<max_BEAGLE_mass_cut)&(start_mass>min_BEAGLE_mass_cut)&(end_mass<max_BEAGLE_mass_cut)&(end_mass>min_BEAGLE_mass_cut)&(start_sfr<max_BEAGLE_sfr_cut)&(start_sfr>min_BEAGLE_sfr_cut)&(end_sfr<max_BEAGLE_sfr_cut)&(end_sfr>min_BEAGLE_sfr_cut)]


start_mass = start_mass[tuple(idx_cut)]
start_sfr = start_sfr[tuple(idx_cut)]
end_mass = end_mass[tuple(idx_cut)]
end_sfr =end_sfr[tuple(idx_cut)]


# =============================================================================
# magnification, 0 is mag applied, 1 is pre-mag correction
# =============================================================================

mag = data['mag_AD'][tuple(idx_cut)]*0 
start_mass = start_mass+mag
start_sfr = start_sfr+mag
end_mass = end_mass+mag
end_sfr = end_sfr+mag



# =============================================================================
# PLOT 
# =============================================================================


x = np.linspace(min_BEAGLE_mass_cut, max_BEAGLE_mass_cut)

plt.figure(figsize=(10, 10))
plt.title('${} < z < {}$'.format(zLow, zHigh))
plt.xlim(min_BEAGLE_mass_cut, max_BEAGLE_mass_cut)
plt.ylim(min_BEAGLE_sfr_cut, max_BEAGLE_sfr_cut)

# START
plt.scatter(start_mass, start_sfr, marker='x')
start_fit = np.polyfit(start_mass, start_sfr, 1)
plt.plot(x, start_fit[0]*x + start_fit[1], color='#1f77b4', label='x start', linewidth=1)

# END
plt.scatter(end_mass, end_sfr, marker='+', s=100)
end_fit = np.polyfit(end_mass, end_sfr, 1)
plt.plot(x, end_fit[0]*x + end_fit[1], color='#ff7f0e', label='+ end', linewidth=1)

# CONNECT
for i in range(len(end_mass)):
    plt.plot((start_mass[i], end_mass[i]), (start_sfr[i], end_sfr[i]), color='k', alpha=0.1)
    
# SANTINI
plt.plot(x, beta_san[san_bin]*x + alpha_san[san_bin], color='k', label='Santini+17', linewidth=1)

plt.legend()
plt.show()













