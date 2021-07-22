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
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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

subfolder = 'delayed_uniform_logtau'
#subfolder = 'danger_delayed_uniform_logtau'
#subfolder = 'danger_constant'
#subfolder = 'safety'

if subfolder == 'danger_constant':
    mass_sfr_option = '_mStar'
    sfr_option = 'log_sfr_arr'
else:
    mass_sfr_option = '_mStar_delayed'
    sfr_option = 'log_sfr_instant_arr'
    
    
    
# =============================================================================
# SCENARIOS
# =============================================================================

# 1 -> 1,2,3,4,5,6 # SANTINI selection
# 2 -> 7 # BEAGLE selection
# 3 -> 7 # addition of chi2 and sfr cut
# 4 -> 7 # trying to replicate my initial results, requires 2 redshift filters, and original mass cuts, and ABmag H cut (instead of flux)
# 5 -> 7 # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut
# 5 -> 7 # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut

scenarioA = 23
scenarioB = 7

# SANTINI BIN (affects redshift + upper/lower mass cuts + santini MS plot)
#san_bins = [0, 1, 2, 3, 4] # 1.65, 2.5, 3.5, 4.5, 5.5
san_bins = [3]

# =============================================================================
# OPTIONS
# =============================================================================

# PLOTS
plot_BEAGLE_vs_AD_redshift = 1
plot_input_to_2sigma = 1
plot_input_to_2sigma_heatplot = 1 # HAS TO USE BEAGLE VALUES OF STELLAR MASS AND INSTANT SFR
plot_MS = 1
plot_MS_zoom = 0 # only includes MCMC adjusted santini
plot_bootstrap = 0
print_bootstrap = 0


# =============================================================================
# Santini
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

alpha_san = B_san - 9.7*A_san
beta_san = A_san
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_err_san = A_err_san

# =============================================================================
# Santini+17 Original values - obtained by eye - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san0 = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san0 = np.array((1.05, 1.1, 0.9, 0.75, 0.55))
A_err_san0 = np.array((0.03, 0.03, 0.04, 0.05, 0.18))
B_san0 = np.array((1.0, 1.15, 1.25, 1.2, 1.5))
B_err_san0 = np.array((0.05, 0.03, 0.03, 0.06, 0.12))

alpha_san0 = B_san0 - 9.7*A_san0
alpha_err_san0 = (B_err_san0**2 + (9.7*A_err_san0)**2) ** 0.5
beta_san0 = A_san0
beta_err_san0 = A_err_san0


# =============================================================================
# delayed uniform logtau - Hogg on scenario 7
# =============================================================================
alpha_7_k0_d_Hogg = np.array([-7.697, -7.104, -6.291, -6.362, -11.383])
beta_7_k0_d_Hogg = np.array([0.91, 0.825, 0.762, 0.84, 1.463])


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
    
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)
    mag_GMM = np.array([np.log10(AD['MAGNIF'])]*3).transpose() # this is genius
    #print(AD.dtype.names)
    
#    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
    sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/' # real Santini values

#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/saved_astrodeep_pickle/astrodeep_pickle.p','r'))
#    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/{}_pickle.p'.format(subfolder, subfolder),'r'))
    
    # for delayed original only, with new chi2 values added:
    astrodeep_pickle = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p','r'))

 

#    MY SANTINI VALUES
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy' # emma technique I think
    
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method
    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method

#    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_temp3.npy'
    sfr_SAN_beta_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files_santini/sfr_santini_beta_temp3.npy'    
    
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
                    'amp_GMM_3d':           astrodeep_pickle['amp_GMM_3d'],
                    
                    'id_SANTINI':           np.load(sbf+'id_SANTINI.npy', allow_pickle=True).astype(float),
                    'mass_SANTINI':         np.load(sbf+'mass_SANTINI.npy', allow_pickle=True).astype(float),
                    'sfr_SANTINI':          np.load(sbf+'sfr_SANTINI.npy', allow_pickle=True).astype(float),
                    'redshift_SANTINI':     np.load(sbf+'redshift_SANTINI.npy', allow_pickle=True).astype(float),
                    'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float)) # -103 -> nan

                    }

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/data.p'.format(subfolder),'w')) # all redshifts   
    # =============================================================================
    # this was for scenario 6 mass cuts - ended up 0.42 less than santini cuts
    # =============================================================================
#    plt.scatter(data['mass_SANTINI'][data['mass_SANTINI']>-20], data['mass_BEAGLE_stellar'][data['mass_SANTINI']>-20], alpha=0.05)
#    plt.xlabel('SANTINI')
#    plt.ylabel('BEAGLE')
#    plt.plot((5,12),(5,12), color='k')
#    plt.plot((5,12),(5-0.42273974795608815,12-0.42273974795608815), color='r')    
#    plt.show()
#    
#    test = data['mass_SANTINI'][data['mass_SANTINI']>-20] - data['mass_BEAGLE_stellar'][data['mass_SANTINI']>-20]
#    test = test[~np.isnan(test)]
#    print(np.mean(test))

    #%%

    # =============================================================================
    # mass limits + redshift bins
    # =============================================================================
    massLow = [8.3, 8.5, 8.8, 8.8, 8.8]
    massLow = massLow[san_bin]
    
    massHigh = [10.2, 10.6, 10.8, 10.8, 10.8]
    massHigh = massHigh[san_bin]
    
    zLow = [1.3, 2.0, 3.0, 4.0, 5.0]
    zHigh = [2.0, 3.0, 4.0, 5.0, 6.0]
    zLow = zLow[san_bin]
    zHigh = zHigh[san_bin]
    wLim = (np.array(zHigh) - np.array(zLow)) / 2.0
    zLim = zLow + wLim

    print(zLim)

    #Tomczak et al. 2016
    TmassHigh = 9.244 + (0.753*zLim) - (0.090*(zLim**2)) # all galaxies (used by santini)
#    TmassHigh = 9.458 + (0.865*zLim) - (0.132*(zLim**2)) # star forming
    
    
    # =============================================================================
    # SCENARIO A
    # =============================================================================
    
    if scenarioA == 1: # SANTINI selection
        
        idx1 = (data['mass_SANTINI'] + data['mag_SANTINI'] > massLow) # removes all -103.0s and 1x -99
        idx2 = (abs(data['redshift_SANTINI']-zLim) < wLim)
        idx3 = (data['mass_SANTINI'] < massHigh) 
            
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))

    if scenarioA == 2: # BEAGLE selection
        pass        
    if scenarioA == 3: # addition of chi2 and sfr cut
        pass
    if scenarioA == 4: # trying to replicate my initial results, requires 2 redshift filters, and original mass cuts, and ABmag H cut (instead of flux)
        pass
    if scenarioA == 5: # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut
        pass
    if scenarioA == 6: # scenario 5 with adjusted mass cuts to account for IMF (-0.42)
        pass

    if scenarioA == 7: # scenario 6 with GMM 3d not fitted removed
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42))
        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42))
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  
        
    if scenarioA == 8: # scenario 7 with outliers included as Hogg will deal with these
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42)) # mass completeness offset from Santini
        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42)) # mass upper offset from Santini
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  

    if scenarioA == 9: # scenario 8 with ONLY BEAGLE redshift cut
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42)) # mass completeness offset from Santini
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42)) # mass upper offset from Santini
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_1)
#        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  

    if scenarioA == 10: # scenario 7 with ONLY BEAGLE redshift cut
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42)) # mass completeness offset from Santini
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42)) # mass upper offset from Santini
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_1)
#        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  
        
    if scenarioA == 11: # scenario 8 with |redshift_BEAGLE-redshift_AD|<1 and 0.5<redshift_BEAGLE<8.0 and tomczak upper mass cut and FIXED lower mass cut off 8.8-0.42
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (8.8-0.42)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*data['redshift_BEAGLE']) - (0.090*(data['redshift_BEAGLE']**2)) # all galaxies (used by santini)        
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh-0.42)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  

    if scenarioA == 12: # scenario 11 with chi2 cut == 20
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (8.8-0.42)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*data['redshift_BEAGLE']) - (0.090*(data['redshift_BEAGLE']**2)) # all galaxies (used by santini)        
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh-0.42)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<20.0) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))    
        
    if scenarioA == 13: # scenario 11 with lower sfr cut at -5 and upper 10
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut 
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (8.8-0.42)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*data['redshift_BEAGLE']) - (0.090*(data['redshift_BEAGLE']**2)) # all galaxies (used by santini)        
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh-0.42)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))  

    if scenarioA == 14: # scenario 12 with lower sfr cut at -5 and upper 10
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (8.8-0.42)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*data['redshift_BEAGLE']) - (0.090*(data['redshift_BEAGLE']**2)) # all galaxies (used by santini)        
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh-0.42)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<20.0) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))    
        
    if scenarioA == 15: # scenario 14 with Emma lower mass limits (ap 0p8) and Tomczak upper mass fixed for z>4, and chi2 BACK to 9.5.
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = 7.492125 + (0.4355*data['redshift_BEAGLE'])
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow-0.42)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh-0.42)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))    

    if scenarioA == 16: # scenario 15 without the -0.42 (no idea if necessary)
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = 7.492125 + (0.4355*data['redshift_BEAGLE'])
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow)) # mass completeness offset from Santini
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1

        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) # mass upper offset from Santini

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))    
                                

    if scenarioA == 17: # scenario 11 with low mass new Emma, high mass tomczak (and fixed for z>4), -5<sfr<10
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.3287:
                MCLmassLow[i] = 8.1
            elif data['redshift_BEAGLE'][i] > 4.1025:
                MCLmassLow[i] = 9.1
            else:
                MCLmassLow[i] = 6.78710161 + 0.56377905*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow)) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = (data['redshift_BEAGLE'] > 4.0)&(data['redshift_BEAGLE'] < 6.0) # BEAGLE redshift        


        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<2000) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_1)
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_3)
        print(sum(idx))
#        idx = np.logical_and(idx,idx6)
        print(sum(idx))
#        idx = np.logical_and(idx,idx7)
        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
        print(sum(idx))        
#        idx = np.logical_and(idx,idx9)
        print(sum(idx))   
        
    if scenarioA == 18: # visual inspection tests, lower mass is 0.1 more accepting, z>5
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.3287:
                MCLmassLow[i] = 8.1
            elif data['redshift_BEAGLE'][i] > 4.1025:
                MCLmassLow[i] = 9.1
            else:
                MCLmassLow[i] = 6.78710161 + 0.56377905*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) - 0.1) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
#        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
#        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = np.logical_or((data['redshift_BEAGLE'] > 5.0), (data['redshift_AD'] > 5.0)) # BEAGLE redshift        
#

        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
#        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

#        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<2000) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_1)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
        idx = np.logical_and(idx,idx5_3)
        print(sum(idx))
#        idx = np.logical_and(idx,idx6)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx7)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))   



    if scenarioA == 19: # visual inspection tests, lower mass is 0.1 more accepting, z>5
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.3287:
                MCLmassLow[i] = 8.1
            elif data['redshift_BEAGLE'][i] > 4.1025:
                MCLmassLow[i] = 9.1
            else:
                MCLmassLow[i] = 6.78710161 + 0.56377905*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) - 0.1) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
#        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
#        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = ~np.logical_or((data['redshift_BEAGLE'] > 5.0), (data['redshift_AD'] > 5.0)) # BEAGLE redshift        
        idx5_4 = np.logical_or((data['redshift_BEAGLE'] > 3.5), (data['redshift_AD'] > 3.5)) # BEAGLE redshif
#

        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
#        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

#        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<2000) # chi-squared
#        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_1)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
        idx = np.logical_and(idx,idx5_3)
        print(sum(idx))
        idx = np.logical_and(idx,idx5_4)
        print(sum(idx))
#        idx = np.logical_and(idx,idx6)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx7)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx8)
#        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        print(sum(idx))   


    if scenarioA == 20: 
    
    # scenario 11 with low mass new Emma, high mass tomczak (and fixed for z>4), -5<sfr<10
    # lower mass is 0.1 more accepting, if z>3.5 then only visual inspection objects included
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.3287:
                MCLmassLow[i] = 8.1
            elif data['redshift_BEAGLE'][i] > 4.1025:
                MCLmassLow[i] = 9.1
            else:
                MCLmassLow[i] = 6.78710161 + 0.56377905*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (MCLmassLow) - 0.1) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((1.3+8.0)/2.0)) < (((1.3+8.0)/2.0) - 1.3)) # 0.5<redshift_BEAGLE<8.0
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
        
        idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
        idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
        
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)
    
        idx5_zgt3p5 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
            if idx5_zgt3p5_temp:
                idx5_zgt3p5[i] = True        
        idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
        
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

#        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<2000) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], s=5, label='AD, {}'.format(len(data['id_AD'])))
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='clusters and relflag, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='H $<$ 27.5, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='lower mass, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx5_z)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='redshift, {}'.format(sum(idx)))
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_3)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_4)
#        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='upper mass, {}'.format(sum(idx)))
        print(sum(idx))
#        idx = np.logical_and(idx,idx7)
#        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='sfr cut, {}'.format(sum(idx)))
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='GMM 3d cut, {}'.format(sum(idx)))
        print(sum(idx))   
        
        
        plt.xlim(0, 10)
        plt.ylim(0, 10)
        plt.legend()
        plt.show()

        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlim(0, 10)
        plt.ylim(0, 10)
        plt.legend()
        plt.show()

    if scenarioA == 21: 
    
    # scenario 11 with low mass new Emma, high mass tomczak (and fixed for z>4), -5<sfr<10
    # lower mass is 0.1 more accepting, if z>3.5 then only visual inspection objects included
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] < 2.3287:
                MCLmassLow[i] = 8.1
            elif data['redshift_BEAGLE'][i] > 4.1025:
                MCLmassLow[i] = 9.1
            else:
                MCLmassLow[i] = 6.78710161 + 0.56377905*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > MCLmassLow) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((0.5+6.5)/2.0)) < (((0.5+6.5)/2.0) - 0.5)) # 0.5<redshift_BEAGLE<8.0 (only affects up to visual inspection of z=3.5)
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
        
        idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
        idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
        
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)
    
        idx5_zgt3p5 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
            if idx5_zgt3p5_temp:
                idx5_zgt3p5[i] = True        
                
        idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data['redshift_BEAGLE'] < 6.5)

        idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
        
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<30) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103
              
        
        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], s=5, label='AD, {}'.format(len(data['id_AD'])))
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='clusters and relflag, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='H $<$ 27.5, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='lower mass, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx5_z)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='redshift, {}'.format(sum(idx)))
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_3)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_4)
#        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='upper mass, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='chi2, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='sfr cut, {}'.format(sum(idx)))
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='GMM 3d cut, {}'.format(sum(idx)))
        print(sum(idx))   
        
        
        plt.xlim(0, 10)
        plt.ylim(0, 10)
        plt.legend()
        plt.show()

        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlim(0, 10)
        plt.ylim(0, 10)
#        plt.legend()
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('redshift')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['redshift_BEAGLE'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('redshift')
        plt.ylabel('sfr')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('sfr')
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['min_chi2_BEAGLE'][idx], data['new_min_chi2_BEAGLE'][idx], s=5)
        plt.xlabel('original chi2')
        plt.ylabel('new chi2')
        plt.show()
        

    if scenarioA == 22: 
    
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] <2.1789654:
                MCLmassLow[i] = 8.0
            elif data['redshift_BEAGLE'][i] > 4.195:
                MCLmassLow[i] = 9.0
            else:
                MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > MCLmassLow) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((0.5+6.5)/2.0)) < (((0.5+6.5)/2.0) - 0.5)) # 0.5<redshift_BEAGLE<8.0 (only affects up to visual inspection of z=3.5)
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
        
        idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
        idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
        
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)
    
        idx5_zgt3p5 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
            if idx5_zgt3p5_temp:
                idx5_zgt3p5[i] = True        
                
        idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data['redshift_BEAGLE'] < 6.5)

        idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
        
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

        idx7 = (data['new_min_chi2_BEAGLE']>0) & (data['new_min_chi2_BEAGLE']<13.28) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103

        
        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], s=5, label='AD, {}'.format(len(data['id_AD'])))
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='clusters and relflag, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='H $<$ 27.5, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='lower mass, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx5_z)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='redshift, {}'.format(sum(idx)))
        print(sum(idx))
#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_3)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_4)
#        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='upper mass, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='chi2, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='sfr cut, {}'.format(sum(idx)))
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='GMM 3d cut, {}'.format(sum(idx)))
        print(sum(idx))   
        
        
        plt.xlim(0, 10)
        plt.ylim(0, 10)
        plt.legend()
        plt.show()

        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlim(0, 10)
        plt.ylim(0, 10)
#        plt.legend()
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('redshift')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['redshift_BEAGLE'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('redshift')
        plt.ylabel('sfr')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('sfr')
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['min_chi2_BEAGLE'][idx], data['new_min_chi2_BEAGLE'][idx], s=5)
        plt.xlabel('original chi2')
        plt.ylabel('new chi2')
        plt.show()
        

    if scenarioA == 23: 
    
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0) # relflag
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

        idx3_IRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
        idx3_IRAC = np.logical_and(idx3_IRAC, data['redshift_BEAGLE']>4.0)
        idx3_IRAC = ~idx3_IRAC
        
        MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
        for i in range(len(MCLmassLow)):
            if data['redshift_BEAGLE'][i] <2.1789654:
                MCLmassLow[i] = 8.0
            elif data['redshift_BEAGLE'][i] > 4.195:
                MCLmassLow[i] = 9.0
            else:
                MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]

        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > MCLmassLow) 
#        print(MCLmassLow, MCLmassLow-0.1)
#        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > -10.0) 
        
#        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim) # astrodeep redshift
#        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim) # BEAGLE redshift
        idx5_1 = (abs(data['redshift_BEAGLE']-((0.5+6.5)/2.0)) < (((0.5+6.5)/2.0) - 0.5)) # 0.5<redshift_BEAGLE<6.5 (only affects up to visual inspection of z=3.5)
        idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
        idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
        
        idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
        idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
        
        vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)
    
        idx5_zgt3p5 = np.full(len(data['id_AD']), False)
        for i in range(len(data['id_AD'])):
            idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
            if idx5_zgt3p5_temp:
                idx5_zgt3p5[i] = True        
                
        idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data['redshift_BEAGLE'] < 6.5)

        idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
        
        TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
        idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh)) 

        idx7 = (data['new_min_chi2_BEAGLE']>0) & (data['new_min_chi2_BEAGLE']<13.28) # chi-squared
        idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
        
        idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103

        
        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], s=5, label='AD, {}'.format(len(data['id_AD'])))
        
        print(sum(idx1))
        
        idx = np.logical_and(idx1,idx2)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='clusters and relflag, {}'.format(sum(idx)))
        print(sum(idx))
        
        idx = np.logical_and(idx,idx3)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='H $<$ 27.5, {}'.format(sum(idx)))
        print(sum(idx))

        idx = np.logical_and(idx,idx4)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='lower mass, {}'.format(sum(idx)))
        print(sum(idx))

        idx = np.logical_and(idx,idx6)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='upper mass, {}'.format(sum(idx)))
        print(sum(idx))




        idx = np.logical_and(idx,idx5_z)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='low z similar AD and BEAGLE, high z visual, {}'.format(sum(idx)))
        print(sum(idx))
  
#        idx = np.logical_and(idx,idx5_zlt3p5)
#        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='high redshift visual inspection, {}'.format(sum(idx)))
#        print(sum(idx))
#        
#        idx = np.logical_and(idx,idx5_zgt3p5)
#        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='high redshift visual inspection, {}'.format(sum(idx)))
#        print(sum(idx))

        idx = np.logical_and(idx,idx3_IRAC)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='remove z$>$4 if no IRAC, {}'.format(sum(idx)))
        print(sum(idx))

#        idx = np.logical_and(idx,idx5_2)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_3)
#        print(sum(idx))
#        idx = np.logical_and(idx,idx5_4)
#        print(sum(idx))

        idx = np.logical_and(idx,idx7)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='chi2, {}'.format(sum(idx)))
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='sfr cut, {}'.format(sum(idx)))
        print(sum(idx))        
        idx = np.logical_and(idx,idx9)
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5, label='GMM 3d cut, {}'.format(sum(idx)))
        print(sum(idx))   
        
        
        plt.xlim(0, 10)
        plt.ylim(0, 10)
        plt.legend()
        plt.show()

        plt.figure(figsize=(10, 10))
        plt.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlim(0, 10)
        plt.ylim(0, 10)
#        plt.legend()
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['redshift_BEAGLE'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('redshift')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['redshift_BEAGLE'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('redshift')
        plt.ylabel('sfr')
        plt.show()        
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['mass_BEAGLE_stellar'][idx], data['sfr_BEAGLE_instant'][idx], s=5)
        plt.xlabel('mass')
        plt.ylabel('sfr')
        plt.show()
        
        plt.figure(figsize=(5, 5))
        plt.scatter(data['min_chi2_BEAGLE'][idx], data['new_min_chi2_BEAGLE'][idx], s=5)
        plt.xlabel('original chi2')
        plt.ylabel('new chi2')
        plt.show()
        

    
        
    for key in data.keys():
        data[key] = data[key][idx]  

#    print(data['field_AD'], data['id_AD'], data['id_BEAGLE'], data['redshift_BEAGLE'], data['redshift_AD'])

#    plt.scatter(data['mag_AD'], data['mag_SANTINI'], alpha=0.2)
#    plt.plot((-0.3,1.5),(-0.3,1.5), color='k')
#    plt.title('Scenario 23')
#    plt.xlabel('AD mag')
#    plt.ylabel('SANTINI mag')
#    plt.show()
    
    # =============================================================================
    # plot AD vs BEAGLE redshift
    # =============================================================================
    if plot_BEAGLE_vs_AD_redshift == 1:
        low = 0
        high = 10
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], alpha=0.1)
        plt.xlim(low,high)
        plt.ylim(low,high)
        plt.xlabel('zAD')
        plt.ylabel('zBEAGLE')
        plt.plot((low,high),(low,high),color='k')
        plt.show()
    
    # =============================================================================
    # Create INPUT for KELLY
    # =============================================================================

#     SCENARIO 
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_data_z{}.p'.format(subfolder, scenarioA, str(zLow).replace('.','p')),'w'))
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_data_z0p5.p'.format(subfolder, scenarioA),'w')) # all redshifts
   
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_subset_zgt5.p'.format(subfolder, scenarioA),'w')) # all redshifts
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/{}/scenario_{}_subset_zgt3p5_lt5.p'.format(subfolder, scenarioA),'w')) # all redshifts
    
    # =============================================================================
    # SCENARIO B
    # =============================================================================

    if scenarioB == 7:
        temp_mass = data['mass_BEAGLE_stellar']
        temp_sfr = data['sfr_BEAGLE_instant']        


#    mass = temp_mass[(temp_mass>0)&(temp_sfr>-30)] # for scenario 8_7, this and line below removed nothing
#    sfr = temp_sfr[(temp_mass>0)&(temp_sfr>-30)]
    
    mass = temp_mass
    sfr = temp_sfr
    
    # =============================================================================
    # 2 sigma process
    # =============================================================================

    outliers = 1
    
    # values saved as _0 prior to 2sigma clipping
    mass0 = mass
    sfr0 = sfr
    
    fit0 = np.polyfit(mass0, sfr0, 1)
    sfr_residuals0 = sfr0 - (fit0[0]*mass0 + fit0[1])
    sigma0 = np.std(sfr_residuals0)
    
    # =============================================================================
    # some plots
    # =============================================================================
    if plot_input_to_2sigma == 1:
        plt.scatter(mass0, sfr0)
        plt.show()
    print('length of temp_mass', len(temp_mass))
    print('length of 2sigma input', len(mass0))

    # =============================================================================
    # heatplots
    # =============================================================================
    #%%
    if plot_input_to_2sigma_heatplot == 1:
        heatplot_mass = np.array([])
        heatplot_sfr = np.array([])
        for i in range(8):

            IDs_per_field = data['id_BEAGLE'][np.isin(data['field_AD'], i)].astype(int) # IDs per field

            heatplot_sample_IDs = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/sample_outputs/{}{}_BEAGLE_samples.p'.format(subfolder, i, mass_sfr_option),'r'))['ID'].astype(int)  
            
            heatplot_sample_masses = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/sample_outputs/{}{}_BEAGLE_samples.p'.format(subfolder, i, mass_sfr_option),'r'))['log_mStar_arr'][np.isin(heatplot_sample_IDs, IDs_per_field)]
            for h in heatplot_sample_masses:
                heatplot_mass = np.concatenate((heatplot_mass, h))
#            
            heatplot_sample_sfrs = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/{}/sample_outputs/{}{}_BEAGLE_samples.p'.format(subfolder, i, mass_sfr_option),'r'))[sfr_option][np.isin(heatplot_sample_IDs, IDs_per_field)]
            for h in heatplot_sample_sfrs:
                heatplot_sfr = np.concatenate((heatplot_sfr, h))

                
        heatplot_mass = np.array(heatplot_mass)
        heatplot_sfr = np.array(heatplot_sfr)
            
#        plt.scatter(heatplot_mass, heatplot_sfr)
#        plt.show()
            
        plt.figure(figsize=(10, 6))
        xlow = 7
        xhigh = 12
        ylow = -4
        yhigh = 4
#        plt.title('${} < z < {}$'.format(zLow, zHigh), size = 20)
        plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
        plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
        plt.hist2d(heatplot_mass, heatplot_sfr, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
        plt.colorbar()
        cmap = cm.get_cmap('viridis')
        rgba = cmap(0)
        ax = plt.axes()
        ax.set_facecolor(rgba)
        plt.show()
        #%%    
#        num = 1100
#        # subset of points
#        plt.figure(figsize=(10, 6))
#        xlow = 7
#        xhigh = 12
#        ylow = -4
#        yhigh = 4
##        plt.title('${} < z < {}$'.format(zLow, zHigh), size = 20)
#        plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
#        plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
#        plt.hist2d(heatplot_mass[num:num+100], heatplot_sfr[num:num+100], bins=150, norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
#        plt.colorbar()
#        cmap = cm.get_cmap('viridis')
#        rgba = cmap(0)
#        ax = plt.axes()
##        ax.set_facecolor(rgba)
#        plt.show()

        #%%
    # =============================================================================
    # continuing with 2sigma
    # =============================================================================
    while outliers > 0:
        
        fit = np.polyfit(mass, sfr, 1)
        sfr_residuals= sfr - (fit[0]*mass + fit[1])  
        sigma = np.std(sfr_residuals)
        idx = (abs(sfr_residuals)<2.0*sigma)    
        outliers = len(mass) - sum(idx)
        mass = mass[idx]
        sfr = sfr[idx]
        
    fit_clipping = fit

    # =============================================================================
    # FIRST YEAR REPORT
    # =============================================================================
    if plot_MS == 1:
        x = np.linspace(7, 10.5)
        plt.figure(figsize=(6, 6))
#        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel(r'$M_\mathrm{tot}$')
        plt.ylabel(r'$\Psi$')
        plt.scatter(mass0, sfr0, marker='x')
        plt.scatter(mass, sfr, marker='x')
        plt.plot(x, fit0[0]*x + fit0[1], color='#1f77b4', label='Without clipping', linewidth=2)
        plt.plot(x, fit[0]*x + fit[1], color='#d62728', label='With clipping', linewidth=2)
        plt.plot(x, beta_san[san_bin]*x + alpha_san[san_bin], color='#2ca02c', label='Santini+17 MCMC', linewidth=2)
        plt.plot(x, beta_san0[san_bin]*x + alpha_san0[san_bin], color='k', label='Santini+17 Raw', linewidth=2, linestyle=':')
#        plt.plot((6.8,11.9),(-1.4,3.2),linestyle=':')
                 
    #    plt.plot(x, (beta_san[san_bin]-beta_err_san[san_bin])*x + (alpha_san[san_bin]-alpha_err_san[san_bin]), color='#2ca02c', linewidth=2, linestyle=':')
    #    plt.plot(x, (beta_san[san_bin]+beta_err_san[san_bin])*x + (alpha_san[san_bin]+alpha_err_san[san_bin]), color='#2ca02c', linewidth=2, linestyle=':')
                 
        # report results
#        plt.plot(x, 0.84*x - 7.16, color='#9467bd', label='Report', linewidth=2) 
                 
        #plt.title(field.replace('_',''))
        
        
        plt.plot(x, beta_7_k0_d_Hogg[san_bin]*x + alpha_7_k0_d_Hogg[san_bin], color='purple', label='Scenario 7 k0 d Hogg', linewidth=2)
        
        
        
        plt.xlim(6, 12)
        plt.ylim(-3, 5)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
        plt.show()
        print('Without clipping {:1.3g} {:1.3g}'.format(fit0[1],fit0[0]))
        print('With clipping {:1.3g} {:1.3g}'.format(fit[1],fit[0]))
        print('Santini+17 MCMC {:1.3g} {:1.3g}'.format(alpha_san[san_bin],beta_san[san_bin]))
        print('Santini+17 Raw {:1.3g} {:1.3g}'.format(alpha_san0[san_bin],beta_san0[san_bin]))
    

    
    # =============================================================================
    # FIRST YEAR REPORT - ZOOMED
    # =============================================================================
    if plot_MS_zoom == 1:
        x = np.linspace(7, 10.5)
        plt.figure(figsize=(6, 6))
        plt.xlabel(r'$M_\mathrm{tot}$')
        plt.ylabel(r'$\Psi$')
        #plt.scatter(mass0, sfr0, marker='x')
        #plt.scatter(mass, sfr, marker='x')
        plt.plot(x, fit0[0]*x + fit0[1], color='#1f77b4', label='Without clipping', linewidth=2)
        plt.plot(x, fit[0]*x + fit[1], color='#d62728', label='With clipping', linewidth=2)
        plt.plot(x, beta_san[san_bin]*x + alpha_san[san_bin], color='#2ca02c', label='Santini+17', linewidth=2)
                 
        plt.plot(x, (beta_san[san_bin]-beta_err_san[san_bin])*x + (alpha_san[san_bin]-alpha_err_san[san_bin]), color='#2ca02c', linewidth=2, linestyle=':')
        plt.plot(x, (beta_san[san_bin]+beta_err_san[san_bin])*x + (alpha_san[san_bin]+alpha_err_san[san_bin]), color='#2ca02c', linewidth=2, linestyle=':')
                 
        
        # report results
        plt.plot(x, 0.84*x - 7.16, color='#9467bd', label='Report', linewidth=2) 
                 
        #plt.title(field.replace('_',''))
        plt.xlim(8.5, 9.0)
        plt.ylim(-0.3, 0.4)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
        plt.show()
    
    # =============================================================================
    # bootstrapping
    # =============================================================================
    
    mass = mass0
    sfr = sfr0
    
    iterations = 20 # 10000 for report
    samples = len(mass)      # len(mass)== 1310, 1000 for report
    
    alpha = []
    beta = []
    difference_8 = [] # difference between best fit and Santini at mass 10^8
    difference_10 = []
    
    for i in range(iterations):
        idx_random = np.random.choice(range(len(mass)), samples, replace=True)
    #    print(np.sort(idx_random))
        mass_bs = mass[idx_random]
        sfr_bs = sfr[idx_random]
        
        outliers = 1
        
        while outliers > 0:
            
            fit = np.polyfit(mass_bs, sfr_bs, 1)
            sfr_residuals= sfr_bs - (fit[0]*mass_bs + fit[1])  
            sigma = np.std(sfr_residuals)
       
            idx = (abs(sfr_residuals)<2.0*sigma)    
            outliers = len(mass_bs) - sum(idx)
            mass_bs = mass_bs[idx]
            sfr_bs = sfr_bs[idx]

            if len(mass_bs) < 5:
                break
            
        if outliers == 0:

            alpha.append(fit[1])
            beta.append(fit[0])
            difference_8.append( (fit[1]+(8.0*fit[0])) - (alpha_san[san_bin]+(8.0*beta_san[san_bin])) )
            difference_10.append( (fit[1]+(10.0*fit[0])) - (alpha_san[san_bin]+(10.0*beta_san[san_bin])) )
    
    
    
    # =============================================================================
    # FIRST YEAR REPORT
    # =============================================================================
    if plot_bootstrap == 1:
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel(r'$\alpha$')
        plt.ylabel('Count')
        y, x, _ = plt.hist(alpha, bins=20)
        
        plt.plot((np.median(alpha), np.median(alpha)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
        plt.plot((np.percentile(alpha, 16), np.percentile(alpha, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
        plt.plot((np.percentile(alpha, 84), np.percentile(alpha, 84)), (0, y.max()), color='#ff7f0e', linewidth=4, linestyle=':')
        plt.plot((alpha_san[san_bin], alpha_san[san_bin]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
        plt.plot((alpha_san[san_bin]-alpha_err_san[san_bin], alpha_san[san_bin]-alpha_err_san[san_bin]), (0, y.max()), label='68', color='#2ca02c', linewidth=4, linestyle=':')
        plt.plot((alpha_san[san_bin]+alpha_err_san[san_bin], alpha_san[san_bin]+alpha_err_san[san_bin]), (0, y.max()), color='#2ca02c', linewidth=4, linestyle=':')
        plt.plot((fit_clipping[1], fit_clipping[1]), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
        plt.plot((-7.16, -7.16), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_alpha_bootstrap.png')
        plt.show()
        
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel(r'$\beta$')
        plt.ylabel('Count')
        y, x, _ = plt.hist(beta, bins=20)
        
        plt.plot((np.median(beta), np.median(beta)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
        plt.plot((np.percentile(beta, 16), np.percentile(beta, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
        plt.plot((np.percentile(beta, 84), np.percentile(beta, 84)), (0, y.max()), linewidth=4, linestyle=':')
        plt.plot((beta_san[san_bin], beta_san[san_bin]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
        plt.plot((beta_san[san_bin]-beta_err_san[san_bin], beta_san[san_bin]-beta_err_san[san_bin]), (0, y.max()), label='68', color='#2ca02c', linewidth=4, linestyle=':')
        plt.plot((beta_san[san_bin]+beta_err_san[san_bin], beta_san[san_bin]+beta_err_san[san_bin]), (0, y.max()), color='#2ca02c', linewidth=4, linestyle=':')
        plt.plot((fit_clipping[0], fit_clipping[0]), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
        plt.plot((0.84, 0.84), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
        plt.show()
        
        
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel('difference8')
        plt.ylabel('Count')
        y, x, _ = plt.hist(difference_8, bins=20)
        
        plt.plot((np.median(difference_8), np.median(difference_8)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
        plt.plot((np.percentile(difference_8, 16), np.percentile(difference_8, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
        plt.plot((np.percentile(difference_8, 84), np.percentile(difference_8, 84)), (0, y.max()), linewidth=4, linestyle=':')
        
        plt.plot(((fit_clipping[1]+(8.0*fit_clipping[0])) - (alpha_san[san_bin]+(8.0*beta_san[san_bin])), (fit_clipping[1]+(8.0*fit_clipping[0])) - (alpha_san[san_bin]+(8.0*beta_san[san_bin]))), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
        plt.plot(((-7.16+(8.0*0.84)) - (alpha_san[san_bin]+(8.0*beta_san[san_bin])), (-7.16+(8.0*0.84)) - (alpha_san[san_bin]+(8.0*beta_san[san_bin]))), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
        plt.show()
        
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel('difference10')
        plt.ylabel('Count')
        y, x, _ = plt.hist(difference_10, bins=20)
        
        plt.plot((np.median(difference_10), np.median(difference_10)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
        plt.plot((np.percentile(difference_10, 16), np.percentile(difference_10, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
        plt.plot((np.percentile(difference_10, 84), np.percentile(difference_10, 84)), (0, y.max()), linewidth=4, linestyle=':')
        
        plt.plot(((fit_clipping[1]+(10.0*fit_clipping[0])) - (alpha_san[san_bin]+(10.0*beta_san[san_bin])), (fit_clipping[1]+(10.0*fit_clipping[0])) - (alpha_san[san_bin]+(10.0*beta_san[san_bin]))), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
        plt.plot(((-7.16+(10.0*0.84)) - (alpha_san[san_bin]+(10.0*beta_san[san_bin])), (-7.16+(10.0*0.84)) - (alpha_san[san_bin]+(10.0*beta_san[san_bin]))), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
        plt.show()
    
    
    # =============================================================================
    # print results
    # =============================================================================
    
    alphas.append(fit_clipping[1])
    betas.append(fit_clipping[0])

    alphas_bs_16.append(np.percentile(alpha, 16))
    alphas_bs_50.append(np.median(alpha))
    alphas_bs_84.append(np.percentile(alpha, 84))
    
    betas_bs_16.append(np.percentile(beta, 16))
    betas_bs_50.append(np.median(beta))
    betas_bs_84.append(np.percentile(beta, 84))
    

if print_bootstrap == 1:
    print('alphas', alphas)
    print('betas', betas)
    
    print('alphas_bs_16', alphas_bs_16)
    print('alphas_bs_50', alphas_bs_50)
    print('alphas_bs_84', alphas_bs_84)
    
    print('betas_bs_16', betas_bs_16)
    print('betas_bs_50', betas_bs_50)
    print('betas_bs_84', betas_bs_84)




