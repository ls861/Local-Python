#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""



import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle


    
AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)

sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'

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
                'metallicity_BEAGLE':   np.load(sbf+'metallicity_BEAGLE.npy')                
    
                }




low=5
high=11
plt.scatter(data['mass_AD_neb'], data['mass_AD'], alpha=0.01, c=data['mag_AD'], vmin=-0.5, vmax=2.0)
plt.xlabel('ADneb')
plt.ylabel('AD')
plt.xlim(low, high)
plt.ylim(low, high)
plt.plot((low,high),(low,high),color='k')
#plt.colorbar()
plt.show()

plt.scatter(data['mass_AD_neb'], data['mass_BEAGLE_tot'], alpha=0.01, c=data['mag_AD'], vmin=-0.5, vmax=2.0)
plt.xlabel('ADneb')
plt.ylabel('BEAGLEtot')
plt.xlim(low, high)
plt.ylim(low, high)
plt.plot((low,high),(low,high),color='k')
#plt.colorbar()
plt.show()

plt.scatter(data['mass_AD_neb'], data['mass_BEAGLE_stellar'], alpha=0.01, c=data['mag_AD'], vmin=-0.5, vmax=2.0)
plt.xlabel('ADneb')
plt.ylabel('BEAGLEstellar')
plt.xlim(low, high)
plt.ylim(low, high)
plt.plot((low,high),(low,high),color='k')
#plt.colorbar()
plt.show()


# =============================================================================
# plsying with magnification
# =============================================================================

plt.scatter(data['mass_AD_neb']+data['mag_AD'], data['mass_BEAGLE_stellar'], alpha=0.01, c=data['mag_AD'], vmin=-0.5, vmax=2.0)
plt.xlabel('ADneb without mag')
plt.ylabel('BEAGLEstellar')
plt.xlim(low, high)
plt.ylim(low, high)
plt.plot((low,high),(low,high),color='k')
#plt.colorbar()
plt.show()

plt.scatter(data['mass_AD_neb']-data['mag_AD'], data['mass_BEAGLE_stellar'], alpha=0.01, c=data['mag_AD'], vmin=-0.5, vmax=2.0)
plt.xlabel('ADneb with mag x2')
plt.ylabel('BEAGLEstellar')
plt.xlim(low, high)
plt.ylim(low, high)
plt.plot((low,high),(low,high),color='k')
#plt.colorbar()
plt.show()




















