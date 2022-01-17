#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 16:53:40 2021

@author: lester
"""



import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import copy
from astropy.table import Table


fname = "IDs.txt"

UID = np.genfromtxt(fname)[1:]

UID_list = []
for i in range(len(UID[:,0])):
    UID_list.append(str(int(UID[:,1][i])) + '_' + str(int(UID[:,0][i])))

UID_list = np.array(UID_list)
print(UID_list)



# =============================================================================
# NOTES
# =============================================================================

'''
NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

NOTE BEAGLE PARAMS have -101 AD OBJECT WAS NOT A BEAGLE INPUT, -102 AD OBJECT WAS NOT FITTED BY BEAGLE, 
sfr_BEAGLE_instant can also have -30 from during creation of instant sfr

THESE BECOME NAN WHEN TAKING LOG: BEAGLE had -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

Objects not included by SANTINI have -103.0 for all params
Objects not fitted by GMM 2d or 3d are also set to -103.0 just for GMM params

NOTE the -30s, -101s and -102s aren't strict as magnification was added to them!

# log(0) -> -inf (mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb)
# lof(-ve) -> nan (mass_BEAGLE_tot and )

'''

# =============================================================================
# SCENARIOS
# =============================================================================

scenarioA = 32

# =============================================================================
# get data
# =============================================================================
AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
mag_GMM = np.array([np.log10(AD['MAGNIF'])]*3).transpose() # this is genius
print(AD.dtype.names)

#    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/' # real Santini values

#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/saved_astrodeep_pickle/astrodeep_pickle.p','r'))
#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/{}_pickle.p'.format(subfolder, subfolder),'r'))

# for new dust beagle only, ignore new chi2 value:
with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/astrodeep_pickle.p', 'rb') as f:
    astrodeep_pickle = pickle.load(f, encoding='latin1')
    
#ORIGINAL
#astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p','r'))

#    MY SANTINI VALUES
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy' # emma technique I think

#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method

#    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_temp3.npy'
sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_beta_temp3.npy'    

# =============================================================================
# make combined file
# =============================================================================
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
                'RA_AD':                AD['RA'],
                'DEC_AD':               AD['DEC'],
                
                'b_CH1_AD':             AD['b_CH1'],
                'b_errCH1_AD':          AD['b_errCH1'],
                'b_CH2_AD':             AD['b_CH2'],
                'b_errCH2_AD':          AD['b_errCH2'],

                'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), 
                'sfr_SAN_beta':         np.load(sfr_SAN_beta_location), 

                'id_BEAGLE':            astrodeep_pickle['id_BEAGLE'],                     
                'mass_BEAGLE_tot':      np.log10(astrodeep_pickle['mass_BEAGLE_tot']) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_stellar':  np.log10(astrodeep_pickle['mass_BEAGLE_stellar']) - np.log10(AD['MAGNIF']), 
                'sfr_BEAGLE_instant':   astrodeep_pickle['sfr_BEAGLE_instant'] - np.log10(AD['MAGNIF']), 
                'redshift_BEAGLE':      astrodeep_pickle['redshift_BEAGLE'], 
                'redshift_BEAGLE_mean':      astrodeep_pickle['redshift_BEAGLE_mean'], 
                'tau_BEAGLE':           astrodeep_pickle['tau_BEAGLE'], 
                'tauv_BEAGLE':          astrodeep_pickle['tauv_BEAGLE'], 
                'msa_BEAGLE':           astrodeep_pickle['msa_BEAGLE'], 
                'metallicity_BEAGLE':   astrodeep_pickle['metallicity_BEAGLE'],
                'min_chi2_BEAGLE':      astrodeep_pickle['min_chi2_BEAGLE'],
                'new_min_chi2_BEAGLE':  astrodeep_pickle['new_min_chi2_BEAGLE'],
                'Ks_BEAGLE_input':      astrodeep_pickle['Ks'],
                'CH1_BEAGLE_input':     astrodeep_pickle['CH1'],
                'CH2_BEAGLE_input':     astrodeep_pickle['CH2'],

                'ch1_beagle_mag_median':     astrodeep_pickle['ch1_beagle_mag_median'],
                'ch1_beagle_mag_lower':     astrodeep_pickle['ch1_beagle_mag_lower'],
                'ch1_beagle_mag_upper':     astrodeep_pickle['ch1_beagle_mag_upper'],
                'ch2_beagle_mag_median':     astrodeep_pickle['ch2_beagle_mag_median'],
                'ch2_beagle_mag_lower':     astrodeep_pickle['ch2_beagle_mag_lower'],
                'ch2_beagle_mag_upper':     astrodeep_pickle['ch2_beagle_mag_upper'],

                'id_GMM_2d':            astrodeep_pickle['id_GMM_2d'],
                'x_GMM_2d':             astrodeep_pickle['x_GMM_2d'] - mag_GMM,
                'y_GMM_2d':             astrodeep_pickle['y_GMM_2d'] - mag_GMM,
                'xsig_GMM_2d':          astrodeep_pickle['xsig_GMM_2d'],
                'ysig_GMM_2d':          astrodeep_pickle['ysig_GMM_2d'],
                'xycov_GMM_2d':         astrodeep_pickle['xycov_GMM_2d'],
                'amp_GMM_2d':           astrodeep_pickle['amp_GMM_2d'],

                'id_GMM_3d':            astrodeep_pickle['id_GMM_3d'],
                'x_GMM_3d':             astrodeep_pickle['x_GMM_3d'] - mag_GMM,
                'y_GMM_3d':             astrodeep_pickle['y_GMM_3d'] - mag_GMM,
                'z_GMM_3d':             astrodeep_pickle['z_GMM_3d'],
                'xsig_GMM_3d':          astrodeep_pickle['xsig_GMM_3d'],
                'ysig_GMM_3d':          astrodeep_pickle['ysig_GMM_3d'],
                'zsig_GMM_3d':          astrodeep_pickle['zsig_GMM_3d'],
                'xycov_GMM_3d':         astrodeep_pickle['xycov_GMM_3d'],
                'xzcov_GMM_3d':         astrodeep_pickle['xzcov_GMM_3d'],
                'yzcov_GMM_3d':         astrodeep_pickle['yzcov_GMM_3d'],
                'amp_GMM_3d':           astrodeep_pickle['amp_GMM_3d']
                
                # 'id_SANTINI':           np.load(sbf+'id_SANTINI.npy', allow_pickle=True).astype(float),
                # 'mass_SANTINI':         np.load(sbf+'mass_SANTINI.npy', allow_pickle=True).astype(float),
                # 'sfr_SANTINI':          np.load(sbf+'sfr_SANTINI.npy', allow_pickle=True).astype(float),
                # 'redshift_SANTINI':     np.load(sbf+'redshift_SANTINI.npy', allow_pickle=True).astype(float),
                # 'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float)) # -103 -> nan

                }

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p','w')) # all redshifts   



if scenarioA == 32:
    full_UID_list = []
    for i in range(len(data['field_AD'])):
        full_UID_list.append(str(int(data['field_AD'][i])) + '_' + str(int(data['id_BEAGLE'][i])))
    full_UID_list = np.array(full_UID_list)

    idx = np.isin(full_UID_list, UID_list)

    data_new = copy.deepcopy(data)
    for key in data_new.keys():
        data_new[key] = data_new[key][idx] 
    
    print('FINAL: {}'.format(len(data_new[key])))

    # with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/pick_a_peak/pick_a_peak.p', 'wb') as f:
    #     pickle.dump(data_new, f)

    outputTable = Table(data_new)
    outputTable.write('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/pick_a_peak/pick_a_peak.fits', overwrite=True)




















