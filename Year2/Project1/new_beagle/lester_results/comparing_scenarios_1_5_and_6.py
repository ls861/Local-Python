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
san_bins = [2]

# =============================================================================
# OPTIONS
# =============================================================================

# PLOTS
plot_BEAGLE_vs_AD_redshift = 0
plot_input_to_2sigma = 1
plot_input_to_2sigma_heatplot = 0 # HAS TO USE BEAGLE VALUES OF STELLAR MASS AND INSTANT SFR
cr           = 1
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
                    'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), # calculated by me
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

                    'id_SANTINI':           np.load(sbf+'id_SANTINI.npy', allow_pickle=True).astype(float), # all calculated by Santini
                    'mass_SANTINI':         np.load(sbf+'mass_SANTINI.npy', allow_pickle=True).astype(float),
                    'sfr_SANTINI':          np.load(sbf+'sfr_SANTINI.npy', allow_pickle=True).astype(float),
                    'redshift_SANTINI':     np.load(sbf+'redshift_SANTINI.npy', allow_pickle=True).astype(float),
                    'mag_SANTINI':          np.log10(np.load(sbf+'mag_SANTINI.npy', allow_pickle=True).astype(float)) # -103 -> nan

                    }

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
    
    # =============================================================================
    # scenario 1 # SANTINI selection
    # =============================================================================
    idx1 = (data['mass_SANTINI'] + data['mag_SANTINI'] > massLow) # removes all -103.0s and 1x -99
    idx2 = (abs(data['redshift_SANTINI']-zLim) < wLim)
    idx3 = (data['mass_SANTINI'] < massHigh) 
        
#    print(sum(idx1))
    idx = np.logical_and(idx1,idx2)
#    print(sum(idx))
    idx = np.logical_and(idx,idx3)
#    print(sum(idx))

    data_1 = {}
        
    for key in data.keys():
        data_1[key] = data[key][idx] 
        
    # =============================================================================
    # scenario 5 # final Kelly choice, SANTINI mass cuts, double redshift filter, H flux cut
    # =============================================================================
    idx1 = (AD['field']%2.0==0.0) # clusters
    idx2 = (AD['RELFLAG']==1.0)
    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
    
    idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > massLow)
    idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
    idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
    idx6 = (data['mass_BEAGLE_stellar'] < massHigh)
        
    idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
    idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
    
#    print(sum(idx1))
    idx = np.logical_and(idx1,idx2)
#    print(sum(idx))
    idx = np.logical_and(idx,idx3)
#    print(sum(idx))
    idx = np.logical_and(idx,idx4)
#    print(sum(idx))
    idx = np.logical_and(idx,idx5_1)
#    print(sum(idx))
    idx = np.logical_and(idx,idx5_2)
#    print(sum(idx))
    idx = np.logical_and(idx,idx6)
#    print(sum(idx))
    idx = np.logical_and(idx,idx7)
#    print(sum(idx))
    idx = np.logical_and(idx,idx8)
#    print(sum(idx))

    data_5 = {}
        
    for key in data.keys():
        data_5[key] = data[key][idx]         
        
    # =============================================================================
    # scenario 6 # scenario 5 with adjusted mass cuts to account for IMF (-0.42)
    # =============================================================================
    idx1 = (AD['field']%2.0==0.0) # clusters
    idx2 = (AD['RELFLAG']==1.0)
    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
    
    idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42))
    idx5_1 = (abs(AD['ZBEST']-zLim) < wLim)
    idx5_2 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
    idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42))
        
    idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
    idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
    
#    print(sum(idx1))
    idx = np.logical_and(idx1,idx2)
#    print(sum(idx))
    idx = np.logical_and(idx,idx3)
#    print(sum(idx))
    idx = np.logical_and(idx,idx4)
#    print(sum(idx))
    idx = np.logical_and(idx,idx5_1)
#    print(sum(idx))
    idx = np.logical_and(idx,idx5_2)
#    print(sum(idx))
    idx = np.logical_and(idx,idx6)
#    print(sum(idx))
    idx = np.logical_and(idx,idx7)
#    print(sum(idx))
    idx = np.logical_and(idx,idx8)
#    print(sum(idx))

    data_6 = {}
        
    for key in data.keys():
        data_6[key] = data[key][idx]  

    # =============================================================================
    # 3 scenarios same as scenario 6, but OMITTING just one filter, to see what crosses are left
    # =============================================================================

    # MASS
    idx = np.logical_and(idx1,idx2)
    idx = np.logical_and(idx,idx3)
#    idx = np.logical_and(idx,idx4)
    idx = np.logical_and(idx,idx5_1)
    idx = np.logical_and(idx,idx5_2)
    idx = np.logical_and(idx,idx6)
    idx = np.logical_and(idx,idx7)
    idx = np.logical_and(idx,idx8)

    data_mass = {}
    for key in data.keys():
        data_mass[key] = data[key][idx]  

    # REDSHIFT
    idx = np.logical_and(idx1,idx2)
    idx = np.logical_and(idx,idx3)
    idx = np.logical_and(idx,idx4)
    idx = np.logical_and(idx,idx5_1)
#    idx = np.logical_and(idx,idx5_2)
    idx = np.logical_and(idx,idx6)
    idx = np.logical_and(idx,idx7)
    idx = np.logical_and(idx,idx8)
    
    data_redshift = {}
    for key in data.keys():
        data_redshift[key] = data[key][idx]  

    # CHI2
    idx = np.logical_and(idx1,idx2)
    idx = np.logical_and(idx,idx3)
    idx = np.logical_and(idx,idx4)
    idx = np.logical_and(idx,idx5_1)
    idx = np.logical_and(idx,idx5_2)
    idx = np.logical_and(idx,idx6)
#    idx = np.logical_and(idx,idx7)
    idx = np.logical_and(idx,idx8)
       
    data_chi2 = {}
    for key in data.keys():
        data_chi2[key] = data[key][idx]  
        
        
    # =============================================================================
    # scenario 7 # quick test for emma for z<0.5
    # =============================================================================
#    idx1 = (AD['field']==0.0) # clusters
#    idx2 = (AD['RELFLAG']==1.0)
#    idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5)
#    
#    idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > (massLow-0.42))
#    idx5_1 = (AD['ZBEST'] < 0.5)
#    idx5_2 = (AD['ZBEST'] > 0.1)
#    idx6 = (data['mass_BEAGLE_stellar'] < (massHigh-0.42))
#        
#    idx7 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<9.5)
#    idx8 = (data['sfr_BEAGLE_instant']>-2.0) & (data['sfr_BEAGLE_instant']<10.0)
#    
##    print(sum(idx1))
#    idx = np.logical_and(idx1,idx2)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx3)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx4)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx5_1)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx5_2)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx6)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx7)
##    print(sum(idx))
#    idx = np.logical_and(idx,idx8)
##    print(sum(idx))
#
#    data_7 = {}
#        
#    for key in data.keys():
#        data_7[key] = data[key][idx] 
#
#    print(len(data_7['id_AD']))
#
#    print(data_7['redshift_AD'])
#    print(data_7['mass_AD'])
#    print(data_7['id_AD']) 
#    print(data_7['id_BEAGLE']) 
#    print(data_7['field_AD'])    

    
    # =============================================================================
    # plots
    # =============================================================================


    idx_B = ((data_6['mass_SANTINI']>0)&(data_6['sfr_SANTINI']>-30))
    
#    plt.figure(figsize=(5,5))
#    plt.title('${} < z < {}$'.format(zLow, zHigh))
#    plt.scatter(data_1['mass_SANTINI'], data_1['sfr_SANTINI'], label='scenario 1')
#    plt.scatter(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], marker='x', label='scenario 6')
#    #plt.scatter(data_5['mass_SANTINI'], data_5['sfr_SANTINI'])
#    plt.show()


    
    test_1 = data_1['mass_SANTINI'] * data_1['sfr_SANTINI']
    test_6 = data_6['mass_SANTINI'][idx_B] * data_6['sfr_SANTINI'][idx_B]
    t = np.isin(test_1, test_6)
    s = np.isin(test_6, test_1)

    
    # =============================================================================
    # mass - sfr
    # =============================================================================
    
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5)
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2)
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), marker='x', zorder=0)
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

    # =============================================================================
    # mass
    # =============================================================================

    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['mass_BEAGLE_stellar'][t], label='scenario 1 and 6', zorder=1, s=5)
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['mass_BEAGLE_stellar'][~t], label='scenario 1 only', marker='x', zorder=2)
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['mass_BEAGLE_stellar'][idx_B][~s], label='scenario 6 only', marker='x', zorder=0)
    plt.xlabel('mass SANTINI')
    plt.ylabel('mass BEAGLE')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset')    
    plt.legend()
    plt.show()

    # =============================================================================
    # mass pre mag
    # =============================================================================

    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter((data_1['mass_SANTINI'] + data_1['mag_SANTINI'])[t], (data_1['mass_BEAGLE_stellar'] + data_1['mag_AD'])[t], label='scenario 1 and 6', zorder=1, s=5)
    plt.scatter((data_1['mass_SANTINI'] + data_1['mag_SANTINI'])[~t], (data_1['mass_BEAGLE_stellar'] + data_1['mag_AD'])[~t], label='scenario 1 only', marker='x', zorder=2)
    plt.scatter((data_6['mass_SANTINI'] + data_6['mag_SANTINI'])[idx_B][~s], (data_6['mass_BEAGLE_stellar'] + data_6['mag_AD'])[idx_B][~s], label='scenario 6 only', marker='x', zorder=0)
    plt.xlabel('pre mag mass SANTINI')
    plt.ylabel('pre mag mass BEAGLE')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset')    
    plt.legend()
    plt.show()
    
    # =============================================================================
    # redshift  
    # =============================================================================

    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter(data_1['redshift_SANTINI'][t], data_1['redshift_BEAGLE'][t], label='scenario 1 and 6', zorder=1, s=5)
    plt.scatter(data_1['redshift_SANTINI'][~t], data_1['redshift_BEAGLE'][~t], label='scenario 1 only', marker='x', zorder=2)
    plt.scatter(data_6['redshift_SANTINI'][idx_B][~s], data_6['redshift_BEAGLE'][idx_B][~s], label='scenario 6 only', marker='x', zorder=0)
    plt.xlabel('redshift SANTINI')
    plt.ylabel('redshift BEAGLE')
    plt.xlim(1, 6)
    plt.ylim(1, 6)
    plt.plot((3,4),(3,3), color='k')
    plt.plot((3,4),(4,4), color='k')
    plt.plot((1,6),(1,6), color='k')
    plt.legend()
    plt.show()


    # =============================================================================
    # colour coding by redshift on the mass vs mass plot
    # =============================================================================

#    Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r

    plt.figure(figsize=(9,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))

    from pylab import cm
    cmap = cm.get_cmap('viridis',3)    # 11 discrete colors
#    plt.scatter(data_1['mass_SANTINI'][t], data_1['mass_BEAGLE_stellar'][t], label='scenario 1 and 6', zorder=1, s=5, cmap=cmap, c=data_1['redshift_BEAGLE'][t], vmin=zLow-1, vmax=zHigh+1)
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['mass_BEAGLE_stellar'][~t], label='scenario 1 only', zorder=2, marker='v', cmap=cmap, c=data_1['redshift_BEAGLE'][~t], vmin=zLow-(zHigh-zLow), vmax=zHigh+(zHigh-zLow))
#    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['mass_BEAGLE_stellar'][idx_B][~s], label='scenario 6 only', zorder=0, marker='x', cmap=cmap, c=data_6['redshift_BEAGLE'][idx_B][~s], vmin=zLow-1, vmax=zHigh+1)
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('mass BEAGLE')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset')
    plt.colorbar(label='redshift BEAGLE')
    plt.legend()
    plt.show()

    # =============================================================================
    # colour coding by chi2 on the mass vs mass plot
    # =============================================================================

    plt.figure(figsize=(9,7))
    plt.title('${} < z < {}$'.format(zLow, zHigh))

    cmap = cm.get_cmap('viridis',3)    # 11 discrete colors
#    plt.scatter(data_1['mass_SANTINI'][t], data_1['mass_BEAGLE_stellar'][t], label='scenario 1 and 6', zorder=1, s=5, cmap=cmap, c=data_1['redshift_BEAGLE'][t], vmin=zLow-1, vmax=zHigh+1)
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['mass_BEAGLE_stellar'][~t], label='scenario 1 only', zorder=2, marker='v', cmap=cmap, c=data_1['min_chi2_BEAGLE'][~t], vmin=0, vmax=28.5)
#    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['mass_BEAGLE_stellar'][idx_B][~s], label='scenario 6 only', zorder=0, marker='x', cmap=cmap, c=data_6['redshift_BEAGLE'][idx_B][~s], vmin=zLow-1, vmax=zHigh+1)
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('mass BEAGLE')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset')
    plt.colorbar(label='chi2 BEAGLE')
    plt.legend()
    plt.show()

    # =============================================================================
    # colour coding by chi2 on the mass vs sfr plot
    # =============================================================================

    viridis = cm.get_cmap('viridis')
    
    plt.figure(figsize=(7,5))
    plt.title('${} < z < {}$'.format(zLow, zHigh))

    cmap = cm.get_cmap('viridis',3)    # 11 discrete colors


    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
#    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2)
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), s=5, zorder=0, color='r')


    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only', zorder=2, marker='x', cmap=cmap, c=data_1['min_chi2_BEAGLE'][~t], vmin=0, vmax=28.5)
    
    
    x = np.linspace(6, 11, 1000)
    
    fit = np.polyfit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], 1)
    plt.plot(x, x*fit[0] + fit[1], color='k')

#    fit = np.polyfit(np.concatenate((data_1['mass_SANTINI'][t],data_6['mass_SANTINI'][idx_B][~s])), np.concatenate((data_1['sfr_SANTINI'][t],data_6['sfr_SANTINI'][idx_B][~s])), 1)
#    plt.plot(x, x*fit[0] + fit[1], color='r')    
    
    fit = np.polyfit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], color='r')    

    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    fit = np.polyfit(data_mass['mass_SANTINI'][idx_mass], data_mass['sfr_SANTINI'][idx_mass], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.0)) 
    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    fit = np.polyfit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.5)) 
    
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    fit = np.polyfit(data_chi2['mass_SANTINI'][idx_chi2], data_chi2['sfr_SANTINI'][idx_chi2], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(1.0))   
    
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.colorbar(label='chi2 BEAGLE')
    plt.legend()
    plt.show()
    
    
    # these are ALL CROSSES 
#    data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t]
    
    idx_B = ((data_6['mass_SANTINI']>0)&(data_6['sfr_SANTINI']>-30))   
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))   

    
    test_1 = data_1['mass_SANTINI'] * data_1['sfr_SANTINI']
    test_6 = data_6['mass_SANTINI'][idx_B] * data_6['sfr_SANTINI'][idx_B]
    test_mass = data_mass['mass_SANTINI'][idx_mass] * data_mass['sfr_SANTINI'][idx_mass]
    test_redshift = data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]
    test_chi2 = data_chi2['mass_SANTINI'][idx_chi2] * data_chi2['sfr_SANTINI'][idx_chi2]
    
    t = np.isin(test_1, test_6)
    s = np.isin(test_6, test_1)
    
    

    t_mass = np.isin(test_1, test_mass)
    t_redshift = np.isin(test_1, test_redshift)
    t_chi2 = np.isin(test_1, test_chi2)
   

    s_mass = np.isin(test_mass, test_1) # applied to data6, gives objects in data1 AND data6
    s_redshift = np.isin(test_redshift, test_1) # applied to data6, gives objects in data1 AND data6
    s_chi2 = np.isin(test_chi2, test_1) # applied to data6, gives objects in data1 AND data6
    
    extra_crosses_mass = np.isin(test_1[t_mass], test_1[~t]) # some crosses
    extra_dots_mass = ~np.isin(test_mass[~s_mass], test_6[~s]) # extra red dots

    extra_crosses_redshift = np.isin(test_1[t_redshift], test_1[~t]) # some crosses
    extra_dots_redshift = ~np.isin(test_redshift[~s_redshift], test_6[~s]) # extra red dots
    
    extra_crosses_chi2 = np.isin(test_1[t_chi2], test_1[~t]) # some crosses
    extra_dots_chi2 = ~np.isin(test_chi2[~s_chi2], test_6[~s]) # extra red dots
    
    # =============================================================================
    # colour coding by redshift on the mass vs sfr plot
    # =============================================================================

    plt.figure(figsize=(7,5))
    plt.title('${} < z < {}$'.format(zLow, zHigh))

    cmap = cm.get_cmap('viridis',3)    # 11 discrete colors


    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
#    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2)
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), s=5, zorder=0, color='r')


    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only', zorder=2, marker='x', cmap=cmap, c=data_1['redshift_BEAGLE'][~t], vmin=zLow-(zHigh-zLow), vmax=zHigh+(zHigh-zLow))
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.colorbar(label='redshift BEAGLE')
    plt.legend()
    plt.show()

    # =============================================================================
    # colour coding by beagle mass on the mass vs sfr plot
    # =============================================================================

    plt.figure(figsize=(7,5))
    plt.title('${} < z < {}$'.format(zLow, zHigh))

    cmap = cm.get_cmap('viridis',4)    # 11 discrete colors


    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
#    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2)
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='r', s=5)


    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only', zorder=2, marker='x', cmap=cmap, c=(data_1['mass_BEAGLE_stellar'] + data_1['mag_AD'])[~t], vmin=(massLow-0.42)-2, vmax=(massLow-0.42)+2)
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.colorbar(label='mass BEAGLE')
    plt.legend()
    plt.show()
    



    # =============================================================================
    # final plot showing mass-sfr, with the extra lines from black, red, new red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2, color=viridis(0.5))
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='r', s=5)
    
    x = np.linspace(6, 11, 1000)
    
    fit = np.polyfit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], 1)
    plt.plot(x, x*fit[0] + fit[1], color='k', label='scenario 1 + 6 (black), {}'.format(len(data_1['mass_SANTINI'][t])))

#    fit = np.polyfit(np.concatenate((data_1['mass_SANTINI'][t],data_6['mass_SANTINI'][idx_B][~s])), np.concatenate((data_1['sfr_SANTINI'][t],data_6['sfr_SANTINI'][idx_B][~s])), 1)
#    plt.plot(x, x*fit[0] + fit[1], color='r')    
    
    fit = np.polyfit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], color='r', label='scenario 6 (red + black), {}'.format(len(data_6['mass_SANTINI'][idx_B])))    

    fit = np.polyfit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.5), label='scenario 1 (green + black), {}'.format(len(data_1['mass_SANTINI'])))   
    
    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    fit = np.polyfit(data_mass['mass_SANTINI'][idx_mass], data_mass['sfr_SANTINI'][idx_mass], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.0), label='scenario 6, no mass cut, {}'.format(len(data_mass['mass_SANTINI'][idx_mass]))) 
    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    fit = np.polyfit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.3), label='scenario 6, no redshift cut, {}'.format(len(data_redshift['mass_SANTINI'][idx_redshift]))) 
    
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    fit = np.polyfit(data_chi2['mass_SANTINI'][idx_chi2], data_chi2['sfr_SANTINI'][idx_chi2], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(1.0), label='scenario 6, no chi2 cut, {}'.format(len(data_chi2['mass_SANTINI'][idx_chi2])))   
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

    # =============================================================================
    # final plot showing mass-sfr, with the extra lines from black, red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2, color=viridis(0.5))
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='r', s=5)
    
    x = np.linspace(6, 11, 1000)
    
    fit = np.polyfit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], 1)
    plt.plot(x, x*fit[0] + fit[1], color='k', label='scenario 1 + 6 (black), {}'.format(len(data_1['mass_SANTINI'][t])))

#    fit = np.polyfit(np.concatenate((data_1['mass_SANTINI'][t],data_6['mass_SANTINI'][idx_B][~s])), np.concatenate((data_1['sfr_SANTINI'][t],data_6['sfr_SANTINI'][idx_B][~s])), 1)
#    plt.plot(x, x*fit[0] + fit[1], color='r')    
    
    fit = np.polyfit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], color='r', label='scenario 6 (red + black), {}'.format(len(data_6['mass_SANTINI'][idx_B])))    

    fit = np.polyfit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.5), label='scenario 1 (green + black), {}'.format(len(data_1['mass_SANTINI'])))   
    
    x_mass = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_mass]))
    y_mass = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_mass]))
    
    x_redshift = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_redshift]))
    y_redshift = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_redshift]))

    x_chi2 = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_chi2]))
    y_chi2 = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_chi2]))    
    
    
    idx_mass = ((x_mass>0)&(y_mass>-30))
    fit = np.polyfit(x_mass[idx_mass], y_mass[idx_mass], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.0), label='scenario 6, no mass cut, {}'.format(len(x_mass[idx_mass]))) 

    idx_redshift = ((x_redshift>0)&(y_redshift>-30))
    fit = np.polyfit(x_redshift[idx_redshift], y_redshift[idx_redshift], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(0.3), label='scenario 6, no redshift cut, {}'.format(len(x_redshift[idx_redshift]))) 
    
    idx_chi2 = ((x_chi2>0)&(y_chi2>-30))
    fit = np.polyfit(x_chi2[idx_chi2], y_chi2[idx_chi2], 1)
    plt.plot(x, x*fit[0] + fit[1], color=viridis(1.0), label='scenario 6, no chi2 cut, {}'.format(len(x_chi2[idx_chi2]))) 
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()



    # =============================================================================
    # now need to compare beagle and santini masses and sfrs
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, {} objects'.format(zLow, zHigh, len(data_6['mass_BEAGLE_stellar'][idx_B])))    
    
    plt.scatter(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], marker='x', color='#d62728')
    fit = np.polyfit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], label='santini mass, santini sfr', color='#d62728')
    
    plt.scatter(data_6['mass_BEAGLE_stellar'][idx_B], data_6['sfr_SANTINI'][idx_B], marker='x', color='#9467bd')
    fit = np.polyfit(data_6['mass_BEAGLE_stellar'][idx_B], data_6['sfr_SANTINI'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], label='beagle mass, santini sfr', color='#9467bd')  

    plt.scatter(data_6['mass_SANTINI'][idx_B], data_6['sfr_BEAGLE_instant'][idx_B], marker='x', color='#8c564b')
    fit = np.polyfit(data_6['mass_SANTINI'][idx_B], data_6['sfr_BEAGLE_instant'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], label='santini mass, beagle sfr', color='#8c564b')  
    
    plt.scatter(data_6['mass_BEAGLE_stellar'][idx_B], data_6['sfr_BEAGLE_instant'][idx_B], marker='x', color='#e377c2')
    fit = np.polyfit(data_6['mass_BEAGLE_stellar'][idx_B], data_6['sfr_BEAGLE_instant'][idx_B], 1)
    plt.plot(x, x*fit[0] + fit[1], label='beagle mass, beagle sfr', color='#e377c2')    
    
    plt.xlabel('mass')
    plt.ylabel('sfr')
    plt.legend()    
    plt.show()    
    

#    print(idx_B)
    #%%
    # =============================================================================
    # above with inclusion of 2 sigma cuts
    # =============================================================================
    def two_sigma_fit(mass_2s, sfr_2s):
        outliers_2s = 1
        while outliers_2s > 0: 
            fit_2s = np.polyfit(mass_2s, sfr_2s, 1)
            sfr_residuals_2s = sfr_2s - (fit_2s[0]*mass_2s + fit_2s[1])  
            sigma_2s = np.std(sfr_residuals_2s)
            idx_2s = (abs(sfr_residuals_2s)<2.0*sigma_2s)    
            outliers_2s = len(mass_2s) - sum(idx_2s)
            mass_2s = mass_2s[idx_2s]
            sfr_2s = sfr_2s[idx_2s]
            count = len(mass_2s)
            
        return fit_2s, count, mass_2s, sfr_2s
    
    
    
    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, {} objects, 2 sigma'.format(zLow, zHigh, len(data_6['mass_BEAGLE_stellar'][idx_B])))    

    mass_2s = data_1['mass_SANTINI']
    sfr_2s = data_1['sfr_SANTINI']
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(mass_2s,sfr_2s)
    plt.scatter(mass_2s, sfr_2s, marker='+', color='#2ca02c', s=80)
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='santini objects, mass and sfr, {}, {}'.format(count, len(data_1['mass_SANTINI'])), color='#2ca02c')          

    mass_2s = data_6['mass_SANTINI'][idx_B]
    sfr_2s = data_6['sfr_SANTINI'][idx_B]
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(mass_2s,sfr_2s)
    plt.scatter(mass_2s, sfr_2s, marker='x', color='#d62728')
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='santini mass, santini sfr, {}'.format(count), color='#d62728')

    mass_2s = data_6['mass_BEAGLE_stellar'][idx_B]
    sfr_2s = data_6['sfr_SANTINI'][idx_B]
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(mass_2s,sfr_2s)
#    plt.scatter(mass_2s, sfr_2s, marker='x', color='#9467bd')
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='beagle mass, santini sfr, {}'.format(count), color='#9467bd', linestyle='dashed')
             
    mass_2s = data_6['mass_SANTINI'][idx_B]
    sfr_2s = data_6['sfr_BEAGLE_instant'][idx_B]
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(mass_2s,sfr_2s)
#    plt.scatter(mass_2s, sfr_2s, marker='x', color='#8c564b')
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='santini mass, beagle sfr, {}'.format(count), color='#8c564b', linestyle='dashed')
             
    mass_2s = data_6['mass_BEAGLE_stellar'][idx_B]
    sfr_2s = data_6['sfr_BEAGLE_instant'][idx_B]
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(mass_2s,sfr_2s)
    plt.scatter(mass_2s, sfr_2s, marker='x', color='#e377c2')
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='beagle mass, beagle sfr, {}'.format(count), color='#e377c2')

 
             
    plt.xlabel('mass')
    plt.ylabel('sfr')
    plt.legend()    
    plt.show()    


    # =============================================================================
    # 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red, new red + some crosses
    # =============================================================================


#%%
    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, 2 sigma'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2, color='#2ca02c')
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)
    
    x = np.linspace(6, 11, 1000)
    
    
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 + 6 (black), {}'.format(count), color='k')

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
    
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black), {}'.format(count), color='#2ca02c')    

    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))    
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_mass['mass_SANTINI'][idx_mass], data_mass['sfr_SANTINI'][idx_mass])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no mass cut, {}'.format(count), color=viridis(0.0))    
    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no redshift cut, {}'.format(count), color=viridis(0.3))
     
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_chi2['mass_SANTINI'][idx_chi2], data_chi2['sfr_SANTINI'][idx_chi2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no chi2 cut, {}'.format(count), color=viridis(1.0))
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

#%%


    # =============================================================================
    # 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, 2 sigma'.format(zLow, zHigh))
    
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=1, s=5, color='k')
    plt.scatter(data_1['mass_SANTINI'][~t], data_1['sfr_SANTINI'][~t], label='scenario 1 only, {}'.format(len(data_1['mass_SANTINI'][~t])), marker='x', zorder=2, color='#2ca02c')
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)
    
    x = np.linspace(6, 11, 1000)
    
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color='k', label='scenario 1 + 6 (black), {}'.format(count))
    
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color='#d62728', label='scenario 6 (red + black), {}'.format(count))    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color='#2ca02c', label='scenario 1 (green + black), {}'.format(count))   
    
    x_mass = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_mass]))
    y_mass = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_mass]))
    
    x_redshift = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_redshift]))
    y_redshift = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_redshift]))

    x_chi2 = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_chi2]))
    y_chi2 = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_chi2]))    
    
    
    idx_mass = ((x_mass>0)&(y_mass>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_mass[idx_mass], y_mass[idx_mass])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color=viridis(0.0), label='scenario 6, no mass cut, {}'.format(count)) 

    idx_redshift = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift], y_redshift[idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color=viridis(0.3), label='scenario 6, no redshift cut, {}'.format(count)) 
    
    idx_chi2 = ((x_chi2>0)&(y_chi2>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_chi2[idx_chi2], y_chi2[idx_chi2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], color=viridis(1.0), label='scenario 6, no chi2 cut, {}'.format(count))
    
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()





    # =============================================================================
    # NO MASS CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red, new red + some crosses
    # =============================================================================


#%%
    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
    plt.scatter(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['sfr_SANTINI'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                 
    plt.scatter(data_1['mass_SANTINI'][~t_mass], data_1['sfr_SANTINI'][~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_mass][extra_crosses_mass], data_1['sfr_SANTINI'][t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_mass['mass_SANTINI'][idx_mass], data_mass['sfr_SANTINI'][idx_mass])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no mass cut (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

#%%


    # =============================================================================
    # NO MASS CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)   
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
#    plt.scatter(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['sfr_SANTINI'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_mass], data_1['sfr_SANTINI'][~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_mass][extra_crosses_mass], data_1['sfr_SANTINI'][t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      

    x_mass = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_mass]))
    y_mass = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_mass]))
    
    idx_mass = ((x_mass>0)&(y_mass>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_mass[idx_mass], y_mass[idx_mass])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no mass cut (red + black (+blue crosses)), {}'.format(count), color='#1f77b4')
    

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()




#%%



    # =============================================================================
    # NO REDSHIFT CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red, new red + some crosses
    # =============================================================================


#%%
    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
    plt.scatter(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['sfr_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_redshift], data_1['sfr_SANTINI'][~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_redshift]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift], data_1['sfr_SANTINI'][t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no redshift cut (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

#%%


    # =============================================================================
    # NO REDSHIFT CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)   
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
#    plt.scatter(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['sfr_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_redshift], data_1['sfr_SANTINI'][~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_redshift]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift], data_1['sfr_SANTINI'][t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      

    x_redshift = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_redshift]))
    y_redshift = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_redshift]))
    
    idx_redshift = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift], y_redshift[idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no redshift cut (red + black (+blue crosses)), {}'.format(count), color='#1f77b4')
    

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()




#%%
    

    # =============================================================================
    # NO CHI2 CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red, new red + some crosses
    # =============================================================================


#%%
    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO CHI2 CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    
    plt.scatter(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], data_chi2['sfr_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], label='extra dots, {}'.format(len(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_chi2], data_1['sfr_SANTINI'][~t_chi2], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_chi2]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2], data_1['sfr_SANTINI'][t_chi2][extra_crosses_chi2], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2])), marker='x', zorder=0, color='#1f77b4')    
 
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_chi2['mass_SANTINI'][idx_chi2], data_chi2['sfr_SANTINI'][idx_chi2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no chi2 cut (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

#%%


    # =============================================================================
    # NO CHI2 CUT 2 SIGMA final plot showing mass-sfr, with the extra lines from black, red + some crosses
    # =============================================================================

    plt.figure(figsize=(10,10))
    plt.title('${} < z < {}$, NO CHI2 CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)   
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    
#    plt.scatter(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], data_chi2['sfr_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], label='extra dots, {}'.format(len(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2])), zorder=0, color='#1f77b4', s=5)                
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_chi2], data_1['sfr_SANTINI'][~t_chi2], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_chi2]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2], data_1['sfr_SANTINI'][t_chi2][extra_crosses_chi2], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2])), marker='x', zorder=0, color='#1f77b4')    

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      

    x_chi2 = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_chi2]))
    y_chi2 = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_chi2]))
    
    idx_chi2 = ((x_chi2>0)&(y_chi2>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_chi2[idx_chi2], y_chi2[idx_chi2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6, no chi2 cut (red + black (+blue crosses)), {}'.format(count), color='#1f77b4')
    

    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()




    # =============================================================================
    # NO MASS CUT 2 SIGMA final plot showing mass-sfr
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
    # POINTS
    
    plt.scatter(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['sfr_SANTINI'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_mass], data_1['sfr_SANTINI'][~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_mass][extra_crosses_mass], data_1['sfr_SANTINI'][t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    


    # LINES

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    x_mass = np.concatenate((data_1['mass_SANTINI'][t], data_1['mass_SANTINI'][t_mass][extra_crosses_mass]))
    y_mass = np.concatenate((data_1['sfr_SANTINI'][t], data_1['sfr_SANTINI'][t_mass][extra_crosses_mass]))
    idx_mass2 = ((x_mass>0)&(y_mass>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_mass[idx_mass2], y_mass[idx_mass2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 subset (black dots +blue crosses), {}'.format(count), color='#2ca02c', linestyle='dashed')          
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 and 6 (black), {}'.format(count), color='k')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_mass['mass_SANTINI'][idx_mass], data_mass['sfr_SANTINI'][idx_mass])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    x_mass = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_mass]))
    y_mass = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_mass]))
    idx_mass2 = ((x_mass>0)&(y_mass>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_mass[idx_mass2], y_mass[idx_mass2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue crosses)), {}'.format(count), color='#1f77b4', linestyle='dashed')
             
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

    
    # =============================================================================
    # NO REDSHIFT CUT 2 SIGMA final plot showing mass-sfr
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
    # POINTS
    
    plt.scatter(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['sfr_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_redshift], data_1['sfr_SANTINI'][~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_redshift]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift], data_1['sfr_SANTINI'][t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    


    # LINES

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    x_redshift = np.concatenate((data_1['mass_SANTINI'][t], data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift]))
    y_redshift = np.concatenate((data_1['sfr_SANTINI'][t], data_1['sfr_SANTINI'][t_redshift][extra_crosses_redshift]))
    idx_redshift2 = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift2], y_redshift[idx_redshift2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 subset (black dots +blue crosses), {}'.format(count), color='#2ca02c', linestyle='dashed')          
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 and 6 (black), {}'.format(count), color='k')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    x_redshift = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_redshift]))
    y_redshift = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_redshift]))
    idx_redshift2 = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift2], y_redshift[idx_redshift2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue crosses)), {}'.format(count), color='#1f77b4', linestyle='dashed')
             
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()

    # =============================================================================
    # NO REDSHIFT CUT 2 SIGMA final plot showing mass-sfr, BEAGLE MASS AND SFR
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
    # POINTS 
    
    plt.scatter(data_redshift['mass_BEAGLE_stellar'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['sfr_BEAGLE_instant'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['mass_BEAGLE_stellar'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['mass_BEAGLE_stellar'][idx_B][~s], data_6['sfr_BEAGLE_instant'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_BEAGLE_stellar'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_BEAGLE_stellar'][t], data_1['sfr_BEAGLE_instant'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_BEAGLE_stellar'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_BEAGLE_stellar'][~t_redshift], data_1['sfr_BEAGLE_instant'][~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['mass_BEAGLE_stellar'][~t_redshift]), len(data_1['mass_BEAGLE_stellar'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_BEAGLE_stellar'][t_redshift][extra_crosses_redshift], data_1['sfr_BEAGLE_instant'][t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['mass_BEAGLE_stellar'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    


    # LINES

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    x_redshift = np.concatenate((data_1['mass_SANTINI'][t], data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift]))
    y_redshift = np.concatenate((data_1['sfr_SANTINI'][t], data_1['sfr_SANTINI'][t_redshift][extra_crosses_redshift]))
    idx_redshift2 = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift2], y_redshift[idx_redshift2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 subset (black dots +blue crosses), {}'.format(count), color='#2ca02c', linestyle='dashed')          
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 and 6 (black), {}'.format(count), color='k')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_redshift['mass_SANTINI'][idx_redshift], data_redshift['sfr_SANTINI'][idx_redshift])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    x_redshift = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_redshift]))
    y_redshift = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_redshift]))
    idx_redshift2 = ((x_redshift>0)&(y_redshift>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_redshift[idx_redshift2], y_redshift[idx_redshift2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue crosses)), {}'.format(count), color='#1f77b4', linestyle='dashed')
             
    plt.xlabel('mass BEAGLE')
    plt.ylabel('sfr BEAGLE')
    plt.legend()
    plt.show()

    # =============================================================================
    # NO CHI2 CUT 2 SIGMA final plot showing mass-sfr
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO CHI2 CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_chi2 = ((data_chi2['mass_SANTINI']>0)&(data_chi2['sfr_SANTINI']>-30))
    
    # POINTS
    
    plt.scatter(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], data_chi2['sfr_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2], label='extra dots, {}'.format(len(data_chi2['mass_SANTINI'][idx_chi2][~s_chi2][extra_dots_chi2])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['sfr_SANTINI'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_chi2], data_1['sfr_SANTINI'][~t_chi2], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_chi2]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2], data_1['sfr_SANTINI'][t_chi2][extra_crosses_chi2], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2])), marker='x', zorder=0, color='#1f77b4')    


    # LINES

    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'], data_1['sfr_SANTINI'])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 (green + black (+blue crosses)), {}'.format(count), color='#2ca02c')  
             
    x_chi2 = np.concatenate((data_1['mass_SANTINI'][t], data_1['mass_SANTINI'][t_chi2][extra_crosses_chi2]))
    y_chi2 = np.concatenate((data_1['sfr_SANTINI'][t], data_1['sfr_SANTINI'][t_chi2][extra_crosses_chi2]))
    idx_chi22 = ((x_chi2>0)&(y_chi2>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_chi2[idx_chi22], y_chi2[idx_chi22])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 subset (black dots +blue crosses), {}'.format(count), color='#2ca02c', linestyle='dashed')          
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_1['mass_SANTINI'][t], data_1['sfr_SANTINI'][t])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 1 and 6 (black), {}'.format(count), color='k')  
             
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_6['mass_SANTINI'][idx_B], data_6['sfr_SANTINI'][idx_B])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black), {}'.format(count), color='#d62728')
      
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(data_chi2['mass_SANTINI'][idx_chi2], data_chi2['sfr_SANTINI'][idx_chi2])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue dots and crosses)), {}'.format(count), color='#1f77b4')

    x_chi2 = np.concatenate((data_6['mass_SANTINI'][idx_B][~s], data_1['mass_SANTINI'][t_chi2]))
    y_chi2 = np.concatenate((data_6['sfr_SANTINI'][idx_B][~s], data_1['sfr_SANTINI'][t_chi2]))
    idx_chi22 = ((x_chi2>0)&(y_chi2>-30))
    fit_clipping_2s, count, mass_2s, sfr_2s = two_sigma_fit(x_chi2[idx_chi22], y_chi2[idx_chi22])
    plt.plot(x, x*fit_clipping_2s[0] + fit_clipping_2s[1], label='scenario 6 (red + black (+blue crosses)), {}'.format(count), color='#1f77b4', linestyle='dashed')
             
    plt.xlabel('mass SANTINI')
    plt.ylabel('sfr SANTINI')
    plt.legend()
    plt.show()
#%%



    # =============================================================================
    # NO MASS CUT mass-mass
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
    plt.scatter(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['mass_BEAGLE_stellar'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['mass_SANTINI'][idx_B][~s], data_6['mass_BEAGLE_stellar'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['mass_SANTINI'][t], data_1['mass_BEAGLE_stellar'][t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['mass_SANTINI'][~t_mass], data_1['mass_BEAGLE_stellar'][~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['mass_SANTINI'][t_mass][extra_crosses_mass], data_1['mass_BEAGLE_stellar'][t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    

             
    plt.xlabel('mass SANTINI')
    plt.ylabel('mass BEAGLE stellar')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset') 
    plt.legend()
    plt.show()

    # =============================================================================
    # NO MASS CUT PRE MAG mass-mass
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
    plt.scatter((data_mass['mass_SANTINI']+data_mass['mag_SANTINI'])[idx_mass][~s_mass][extra_dots_mass], (data_mass['mass_BEAGLE_stellar']+data_mass['mag_AD'])[idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter((data_6['mass_SANTINI']+data_6['mag_SANTINI'])[idx_B][~s], (data_6['mass_BEAGLE_stellar']+data_6['mag_AD'])[idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[~t_mass], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t_mass][extra_crosses_mass], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    

             
    plt.xlabel('pre mag mass SANTINI')
    plt.ylabel('pre mag mass BEAGLE stellar')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset') 
    plt.legend()
    plt.show()

    # =============================================================================
    # NO MASS CUT PRE MAG mass-mass COLOURED by Z
    # =============================================================================
#%%
    plt.figure(figsize=(9,7))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))

    cmap = cm.get_cmap('hsv')   
    vmin=-1
    vmax=1

    plt.scatter((data_mass['mass_SANTINI']+data_mass['mag_SANTINI'])[idx_mass][~s_mass][extra_dots_mass], (data_mass['mass_BEAGLE_stellar']+data_mass['mag_AD'])[idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['mass_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, s=5, cmap=cmap, vmin=vmin, vmax=vmax, c=(data_mass['redshift_SANTINI']-data_mass['redshift_BEAGLE'])[idx_mass][~s_mass][extra_dots_mass]) 
      
    plt.scatter((data_6['mass_SANTINI']+data_6['mag_SANTINI'])[idx_B][~s], (data_6['mass_BEAGLE_stellar']+data_6['mag_AD'])[idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, s=5, cmap=cmap, vmin=vmin, vmax=vmax, c=(data_6['redshift_SANTINI']-data_6['redshift_BEAGLE'])[idx_B][~s])               
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, cmap=cmap, vmin=vmin, vmax=vmax, c=(data_1['redshift_SANTINI']-data_1['redshift_BEAGLE'])[t])                
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[~t_mass], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_mass]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, cmap=cmap, vmin=vmin, vmax=vmax, c=(data_1['redshift_SANTINI']-data_1['redshift_BEAGLE'])[~t_mass])
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t_mass][extra_crosses_mass], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, cmap=cmap, vmin=vmin, vmax=vmax, c=(data_1['redshift_SANTINI']-data_1['redshift_BEAGLE'])[t_mass][extra_crosses_mass])    

    plt.xlabel('pre mag mass SANTINI')
    plt.ylabel('pre mag mass BEAGLE stellar')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset') 
    plt.colorbar(label='santini z - beagle z')
    plt.legend()
    plt.show()    
    
    # =============================================================================
    # NO MASS CUT z-z
    # =============================================================================

#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO MASS CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_mass = ((data_mass['mass_SANTINI']>0)&(data_mass['sfr_SANTINI']>-30))
    
    plt.scatter(data_mass['redshift_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['redshift_BEAGLE'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['redshift_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['redshift_SANTINI'][idx_B][~s], data_6['redshift_BEAGLE'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['redshift_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['redshift_SANTINI'][t], data_1['redshift_BEAGLE'][t], label='scenario 1 and 6, {}'.format(len(data_1['redshift_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['redshift_SANTINI'][~t_mass], data_1['redshift_BEAGLE'][~t_mass], label='scenario 1 only, {}, {}'.format(len(data_1['redshift_SANTINI'][~t_mass]), len(data_1['redshift_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['redshift_SANTINI'][t_mass][extra_crosses_mass], data_1['redshift_BEAGLE'][t_mass][extra_crosses_mass], label='extra crosses, {}'.format(len(data_1['redshift_SANTINI'][t_mass][extra_crosses_mass])), marker='x', zorder=0, color='#1f77b4')    

    plt.xlim(zLow, zHigh)
    plt.ylim(1, 6)
    plt.plot((3,4),(3,3), color='k')
    plt.plot((3,4),(4,4), color='k')
    plt.plot((1,6),(1,6), color='k')             
    plt.xlabel('redshift SANTINI')
    plt.ylabel('redshift BEAGLE stellar')
    plt.legend()
    plt.show()

    # =============================================================================
    # NO Z CUT PRE MAG mass-mass
    # =============================================================================
#%%
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
    plt.scatter((data_redshift['mass_SANTINI']+data_redshift['mag_SANTINI'])[idx_redshift][~s_redshift][extra_dots_redshift], (data_redshift['mass_BEAGLE_stellar']+data_redshift['mag_AD'])[idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['mass_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter((data_6['mass_SANTINI']+data_6['mag_SANTINI'])[idx_B][~s], (data_6['mass_BEAGLE_stellar']+data_6['mag_AD'])[idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['mass_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t], label='scenario 1 and 6, {}'.format(len(data_1['mass_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[~t_redshift], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['mass_SANTINI'][~t_redshift]), len(data_1['mass_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter((data_1['mass_SANTINI']+data_1['mag_SANTINI'])[t_redshift][extra_crosses_redshift], (data_1['mass_BEAGLE_stellar']+data_1['mag_AD'])[t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['mass_SANTINI'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    

             
    plt.xlabel('pre mag mass SANTINI')
    plt.ylabel('pre mag mass BEAGLE stellar')
    plt.plot((6,12),(6,12), color='k', label='1 to 1')
    plt.plot((6,12),(6-0.4227,12-0.4227), color='r', label='-0.42 offset') 
    plt.legend()
    plt.show()
    
    
    # =============================================================================
    # NO Z CUT z-z
    # =============================================================================

#%%
    
    plt.figure(figsize=(7,7))
    plt.title('${} < z < {}$, NO REDSHIFT CUT, 2 sigma fitted lines'.format(zLow, zHigh))
    x = np.linspace(6, 11, 1000)    
    idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))
    
    plt.scatter(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['redshift_BEAGLE'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)             
    plt.scatter(data_6['redshift_SANTINI'][idx_B][~s], data_6['redshift_BEAGLE'][idx_B][~s], label='scenario 6 only, {}'.format(len(data_6['redshift_SANTINI'][idx_B][~s])), zorder=0, color='#d62728', s=5)               
    plt.scatter(data_1['redshift_SANTINI'][t], data_1['redshift_BEAGLE'][t], label='scenario 1 and 6, {}'.format(len(data_1['redshift_SANTINI'][t])), zorder=0, s=5, color='k')                
    plt.scatter(data_1['redshift_SANTINI'][~t_redshift], data_1['redshift_BEAGLE'][~t_redshift], label='scenario 1 only, {}, {}'.format(len(data_1['redshift_SANTINI'][~t_redshift]), len(data_1['redshift_SANTINI'][~t])), marker='x', zorder=0, color='#2ca02c')
    plt.scatter(data_1['redshift_SANTINI'][t_redshift][extra_crosses_redshift], data_1['redshift_BEAGLE'][t_redshift][extra_crosses_redshift], label='extra crosses, {}'.format(len(data_1['redshift_SANTINI'][t_redshift][extra_crosses_redshift])), marker='x', zorder=0, color='#1f77b4')    

    plt.xlim(zLow, zHigh)
    plt.ylim(1, 6)
    plt.plot((zLow,zHigh),(zLow,zLow), color='k')
    plt.plot((zLow,zHigh),(zHigh,zHigh), color='k')
    plt.plot((1,6),(1,6), color='k')             
    plt.xlabel('redshift SANTINI')
    plt.ylabel('redshift BEAGLE stellar')
    plt.legend()
    plt.show()

'''

print(idx_redshift)
print(len(idx_redshift))
print(sum(idx_redshift))

print(len(s_redshift))

print(len(test_redshift))
print(len())
print(len())
print(len())
print(sum(((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30))))

idx_redshift = ((data_redshift['mass_SANTINI']>0)&(data_redshift['sfr_SANTINI']>-30)) # NEEDED
idx_redshift = ((data_redshift['redshift_SANTINI']>0)&(data_redshift['redshift_BEAGLE']>-30)) # used
        
        
    test_1 = data_1['mass_SANTINI'] * data_1['sfr_SANTINI']
    test_6 = data_6['mass_SANTINI'][idx_B] * data_6['sfr_SANTINI'][idx_B]
    test_mass = data_mass['mass_SANTINI'][idx_mass] * data_mass['sfr_SANTINI'][idx_mass]
    test_redshift = data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]
    test_chi2 = data_chi2['mass_SANTINI'][idx_chi2] * data_chi2['sfr_SANTINI'][idx_chi2]
    
    
    print(len(data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]))
    
    print('FFFFFS', len(test_redshift))
    
    t = np.isin(test_1, test_6)
    s = np.isin(test_6, test_1)
    
    

    t_mass = np.isin(test_1, test_mass)
    t_redshift = np.isin(test_1, test_redshift)
    t_chi2 = np.isin(test_1, test_chi2)
   

    s_mass = np.isin(test_mass, test_1) # applied to data6, gives objects in data1 AND data6
    s_redshift = np.isin(test_redshift, test_1) # applied to data6, gives objects in data1 AND data6
    s_chi2 = np.isin(test_chi2, test_1) # applied to data6, gives objects in data1 AND data6
    
    
    
# =============================================================================
# trying to fix a bug
# =============================================================================

    plt.scatter(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['redshift_BEAGLE'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5) 


print(len(data_redshift['redshift_SANTINI'][idx_redshift]))
print(len(~s_redshift))
print(len(np.isin(test_redshift, test_1)))
print(len(test_redshift))
print(len(test_1))
print(len(data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]))
print(len(test_redshift))
print(len())
print(len())

print(len(data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]))
    
print('FFFFFS', len(test_redshift))
    
'''

'''
print(data_mass['redshift_SANTINI'][idx_mass][~s_mass])
print(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift])
print(data_mass['redshift_BEAGLE'][idx_mass][~s_mass])
print(data_redshift['redshift_BEAGLE'][idx_redshift][~s_redshift])

print(len(data_mass['redshift_BEAGLE'][idx_mass]))
print(len(~s_mass))


    plt.scatter(data_mass['redshift_SANTINI'][idx_mass][~s_mass][extra_dots_mass], data_mass['redshift_BEAGLE'][idx_mass][~s_mass][extra_dots_mass], label='extra dots, {}'.format(len(data_mass['redshift_SANTINI'][idx_mass][~s_mass][extra_dots_mass])), zorder=0, color='#1f77b4', s=5)  
    plt.scatter(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift], data_redshift['redshift_BEAGLE'][idx_redshift][~s_redshift][extra_dots_redshift], label='extra dots, {}'.format(len(data_redshift['redshift_SANTINI'][idx_redshift][~s_redshift][extra_dots_redshift])), zorder=0, color='#1f77b4', s=5)        
                




    test_1 = data_1['mass_SANTINI'] * data_1['sfr_SANTINI']
    test_6 = data_6['mass_SANTINI'][idx_B] * data_6['sfr_SANTINI'][idx_B]
    test_mass = data_mass['mass_SANTINI'][idx_mass] * data_mass['sfr_SANTINI'][idx_mass]
    test_redshift = data_redshift['mass_SANTINI'][idx_redshift] * data_redshift['sfr_SANTINI'][idx_redshift]
    test_chi2 = data_chi2['mass_SANTINI'][idx_chi2] * data_chi2['sfr_SANTINI'][idx_chi2]
    
    t = np.isin(test_1, test_6)
    s = np.isin(test_6, test_1)

    t_mass = np.isin(test_1, test_mass)
    t_redshift = np.isin(test_1, test_redshift)
    t_chi2 = np.isin(test_1, test_chi2)
   
    s_mass = np.isin(test_mass, test_1) # applied to data6, gives objects in data1 AND data6
    s_redshift = np.isin(test_redshift, test_1) # applied to data6, gives objects in data1 AND data6
    s_chi2 = np.isin(test_chi2, test_1) # applied to data6, gives objects in data1 AND data6                
                
print(len(test_1))
print(len(test_mass))
print(len(s_mass))
print(len(data_mass['mass_SANTINI'][idx_mass]))
print(len(data_mass['sfr_SANTINI'][idx_mass]))
print(len(data_mass['mass_SANTINI'][idx_mass] * data_mass['sfr_SANTINI'][idx_mass]))


'''



#%%
'''


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


'''

