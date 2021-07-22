#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# =============================================================================
# # Santini approx 1711 sources 1.3 < z < 6.0
# =============================================================================

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle

# =============================================================================
# san_bin loop
# =============================================================================

alphas = []
alphas_bs_16 = []
alphas_bs_50 = []
alphas_bs_84 = []

betas = []
betas_bs_16 = []
betas_bs_50 = []
betas_bs_84 = []

for san_bin in san_bins:
    
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/astrodeep/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)
    #print(AD.dtype.names)
    
    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_temp3.npy'
   
    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_temp3.npy'

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
                    'sfr_SAN_beta':         np.load(sfr_SAN_beta_location), 
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
                    'amp_GMM':              np.load(sbf+'amp_GMM.npy'),

                    'id_SANTINI':           np.load(sbf+'id_SANTINI.npy', allow_pickle=True).astype(float),
                    'mass_SANTINI':         np.load(sbf+'mass_SANTINI.npy', allow_pickle=True).astype(float),
                    'sfr_SANTINI':          np.load(sbf+'sfr_SANTINI.npy', allow_pickle=True).astype(float),
                    'redshift_SANTINI':     np.load(sbf+'redshift_SANTINI.npy', allow_pickle=True).astype(float),
                    'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float))

                    }
    
print(len(data['mass_SANTINI']))
print(len(data['mass_SANTINI'][data['field_AD']==0]))  
print(len(data['mass_SANTINI'][(data['field_AD']==0)&(data['id_SANTINI']>0)]))  # 594 sources as confirmed by Santini
        
idx_data = (data['field_AD']==0)&(data['id_SANTINI']>0)
        
for key in data.keys():
    data[key] = data[key][idx_data] 

for name in ['id', 'mass', 'sfr', 'redshift', 'mag']:
    plt.title(name)
    plt.scatter(data['{}_AD'.format(name)], data['{}_SANTINI'.format(name)], alpha=0.1)
    plt.xlabel('AD')
    plt.ylabel('SANTINI')
    plt.show()
        



#%%

low = 5
high = 12
plt.scatter(data['mass_AD_neb'], data['mass_SANTINI'], alpha=0.1)
plt.xlim(low,high)
plt.ylim(low,high)
plt.xlabel('massADneb')
plt.ylabel('massSANTINI')
plt.plot((low,high),(low,high),color='k')
plt.show()

low = -5
high = 5
plt.scatter(data['sfr_SAN'], data['sfr_SANTINI'], alpha=0.1)
plt.xlim(low,high)
plt.ylim(low,high)
plt.xlabel('sfrAD')
plt.ylabel('sfrSANTINI')
plt.plot((low,high),(low,high),color='k')
plt.show()

low = 0
high = 10
plt.scatter(data['redshift_AD'], data['redshift_SANTINI'], alpha=0.1)
plt.xlim(low,high)
plt.ylim(low,high)
plt.xlabel('zAD')
plt.ylabel('zSANTINI')
plt.plot((low,high),(low,high),color='k')
plt.show()

low = -0.5
high = 2
plt.scatter(data['mag_AD'], data['mag_SANTINI'], alpha=0.1)
plt.xlim(low,high)
plt.ylim(low,high)
plt.xlabel('magAD')
plt.ylabel('magSANTINI')
plt.plot((low,high),(low,high),color='k')
plt.show()

low = 0
high = 10
plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], alpha=0.1)
plt.xlim(low,high)
plt.ylim(low,high)
plt.xlabel('zAD')
plt.ylabel('zBEAGLE')
plt.plot((low,high),(low,high),color='k')
plt.show()



#print(len(data['redshift_AD']))









#%%



#print(data['id_SANTINI'][idx_santini_0][data['mag_SANTINI'][idx_santini_0]>1.707])
#print(data['id_AD'][idx_santini_0][data['mag_SANTINI'][idx_santini_0]>1.707])
#print(data['mag_SANTINI'][idx_santini_0][data['mag_SANTINI'][idx_santini_0]>1.707])
#print(data['mag_AD'][idx_santini_0][data['mag_SANTINI'][idx_santini_0]>1.707])
#
#
#print(np.log10(9.138113))


        

        