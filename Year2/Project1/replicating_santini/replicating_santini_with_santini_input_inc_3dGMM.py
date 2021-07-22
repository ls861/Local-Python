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

# =============================================================================
# NOTES
# =============================================================================

'''
NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

NOTE BEAGLE PARAMS have -101 AD OBJECT WAS NOT A BEAGLE INPUT, -102 AD OBJECT WAS NOT FITTED BY BEAGLE, 
sfr_BEAGLE_instant can also have -30 from during creation of instant sfr

THESE BECOME NAN WHEN TAKING LOG: BEAGLE had -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

Objects not included by SANTINI have -103.0 for all params

NOTE the -30s, -101s and -102s aren't strict as magnification was added to them!

# log(0) -> -inf (mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb)
# lof(-ve) -> nan (mass_BEAGLE_tot and )

'''

# =============================================================================
# SCENARIOS
# =============================================================================

# 1 -> 1,2,3,4,5,6 # SANTINI selection
# 2 -> 7 # BEAGLE selection
# 3 -> 7 # addition of chi2 and sfr cut
# 4 -> 7 # trying to replicate my initial results, requires 2 redshift filters, and original mass cuts, and ABmag H cut (instead of flux)
# 5 -> 7 # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut
# 5 -> 7 # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut

scenarioA = 6
scenarioB = 7

# SANTINI BIN (affects redshift + upper/lower mass cuts + santini MS plot)
san_bins = [0, 1, 2, 3, 4] # 1.65, 2.5, 3.5, 4.5, 5.5
#san_bins = [0]

# =============================================================================
# OPTIONS
# =============================================================================

# PLOTS
plot_BEAGLE_vs_AD_redshift = 0
plot_input_to_2sigma = 1
plot_input_to_2sigma_heatplot = 0 # HAS TO USE BEAGLE VALUES OF STELLAR MASS AND INSTANT SFR
plot_MS = 1
plot_MS_zoom = 0 # only includes MCMC adjusted santini
plot_bootstrap = 1
print_bootstrap = 1


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
    
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy' # emma technique I think
    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_temp3.npy' # my 1500, central filter, quoted wavelength method
   
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
                    'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float)) # -103 -> nan

                    }
        
   
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
    wLim = (zHigh - zLow) / 2.0
    zLim = zLow + wLim
    
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
        idx_massLim = (abs(data['redshift_SANTINI']-zLim) < wLim)
        massLow_offset = np.mean(((data['mass_SANTINI'][idx_massLim] + data['mag_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim] + data['mag_AD'][idx_massLim]))[np.isfinite((data['mass_SANTINI'][idx_massLim] + data['mag_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim] + data['mag_AD'][idx_massLim]))]) # (ignore nan values)
        massHigh_offset = np.mean(((data['mass_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim]))[np.isfinite((data['mass_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim]))]) # (ignore nan values)
    
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > massLow-massLow_offset)
        idx5 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh-massHigh_offset) 
            
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3) 
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        
    if scenarioA == 3: # addition of chi2 and sfr cut
        idx_massLim = (abs(data['redshift_SANTINI']-zLim) < wLim)
        massLow_offset = np.mean(((data['mass_SANTINI'][idx_massLim] + data['mag_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim] + data['mag_AD'][idx_massLim]))[np.isfinite((data['mass_SANTINI'][idx_massLim] + data['mag_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim] + data['mag_AD'][idx_massLim]))]) # (ignore nan values)
        massHigh_offset = np.mean(((data['mass_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim]))[np.isfinite((data['mass_SANTINI'][idx_massLim]) - (data['mass_BEAGLE_stellar'][idx_massLim]))]) # (ignore nan values)
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > massLow-massLow_offset)
        idx5 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh-massHigh_offset) 
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
        
        print(sum(idx1))
        idx = np.logical_and(idx1,idx2)
        print(sum(idx))
        idx = np.logical_and(idx,idx3)
        print(sum(idx))
        idx = np.logical_and(idx,idx4)
        print(sum(idx))
        idx = np.logical_and(idx,idx5)
        print(sum(idx))
        idx = np.logical_and(idx,idx6)
        print(sum(idx))
        idx = np.logical_and(idx,idx7)
        print(sum(idx))
        idx = np.logical_and(idx,idx8)
        print(sum(idx))

    if scenarioA == 4: # trying to replicate my initial results, requires 2 redshift filters, and original mass cuts, and ABmag H cut (instead of flux)
        massLow_offset = np.mean(((data['mass_AD_neb'] + np.log10(AD['MAGNIF'])) - (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF'])))[np.isfinite((data['mass_AD_neb'] + np.log10(AD['MAGNIF'])) - (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF'])))]) # (ignore nan values)
        massHigh_offset = np.mean((data['mass_AD_neb'] - data['mass_BEAGLE_stellar'])[np.isfinite(data['mass_AD_neb'] - data['mass_BEAGLE_stellar'])]) # (ignore nan values)
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (AD['H160']<27.5)
#        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > massLow-massLow_offset)
        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh-massHigh_offset)
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
        
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

    if scenarioA == 5: # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > massLow)
        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh)
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
        
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

    if scenarioA == 6: # scenario 5 with adjusted mass cuts to account for IMF (-0.42)
        
        idx1 = (AD['field']%2.0==0.0) # clusters
        idx2 = (AD['RELFLAG']==1.0)
        idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
        
        idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42))
        idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
        idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42))
            
        idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
        idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
        
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
        
        
    for key in data.keys():
        data[key] = data[key][idx]  

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

    # SCENARIO 
#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data/scenario_{}_{}_data_z{}.p'.format(scenarioA, scenarioB, str(zLow).replace('.','p')),'w'))
    
    # =============================================================================
    # SCENARIO B
    # =============================================================================

    if scenarioB == 1:
        temp_mass = data['mass_SANTINI']
        temp_sfr = data['sfr_SANTINI']

    if scenarioB == 2:
        temp_mass = data['mass_AD']
        temp_sfr = data['sfr_SANTINI']
        
    if scenarioB == 3:
        temp_mass = data['mass_AD_neb']
        temp_sfr = data['sfr_SANTINI']
        
    if scenarioB == 4:
        temp_mass = data['mass_BEAGLE_stellar']
        temp_sfr = data['sfr_SANTINI']
    
    
    
    if scenarioB == 5:
        temp_mass = data['mass_SANTINI']
        temp_sfr = data['sfr_SAN']    
    
    if scenarioB == 6:
        temp_mass = data['mass_SANTINI']
        temp_sfr = data['sfr_BEAGLE_instant']    
        
        
        
    if scenarioB == 7:
        temp_mass = data['mass_BEAGLE_stellar']
        temp_sfr = data['sfr_BEAGLE_instant']        

    print('length of temp_mass', len(temp_mass))
    mass = temp_mass[(temp_mass>0)&(temp_sfr>-30)]
    sfr = temp_sfr[(temp_mass>0)&(temp_sfr>-30)]
    
    
    
    
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
    print('length of 2sigma input', len(mass0))

    # =============================================================================
    # heatplots
    # =============================================================================
    if plot_input_to_2sigma_heatplot == 1:
        heatplot_mass = np.array([])
        heatplot_sfr = np.array([])
        for i in range(len(mass0)):
            heatplot_samples = pickle.load(open('/Users/lester/Documents/BEAGLE_heatplot_samples/{}_{}_BEAGLE_samples.p'.format(int(data['field_AD'][i]), data['id_BEAGLE'][i]),'r'))
            heatplot_mass = np.concatenate((heatplot_mass, heatplot_samples['mStar']))
            heatplot_sfr = np.concatenate((heatplot_sfr, heatplot_samples['sfr_instant']))
            
        heatplot_mass = np.array(heatplot_mass)
        heatplot_sfr = np.array(heatplot_sfr)
            
        plt.scatter(heatplot_mass, heatplot_sfr)
        plt.show()
        
        plt.hexbin(heatplot_mass, heatplot_sfr, gridsize=[200,200])
        plt.xlim(7, 10)
        plt.ylim(-2, 3)
        plt.tight_layout()
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.show()
            
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
        plt.title('${} < z < {}$'.format(zLow, zHigh))
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
        plt.plot(x, 0.84*x - 7.16, color='#9467bd', label='Report', linewidth=2) 
                 
        #plt.title(field.replace('_',''))
        plt.xlim(7, 10.5)
        plt.ylim(-2, 3)
        plt.legend()
        plt.tight_layout()
        #plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
        plt.show()
    
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




