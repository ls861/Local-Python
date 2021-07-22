#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# =============================================================================
# NOTES
# =============================================================================

'''
NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

NOTE BEAGLE PARAMS have -101 AD OBJECT WAS NOT A BEAGLE INPUT, -102 AD OBJECT WAS NOT FITTED BY BEAGLE, 
sfr_BEAGLE_instant can also have -30 from during creation of instant sfr

THESE BECOME NAN WHEN TAKING LOG: BEAGLE had -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

NOTE the -30s, -101s and -102s aren't strict as magnification was added to them!

# log(0) -> -inf (mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb)
# lof(-ve) -> nan (mass_BEAGLE_tot and )

'''

# =============================================================================
# OPTIONS
# =============================================================================

# H<27.5, relflag==1

# SAME + zLow<redshift_BEAGLE<zHigh
option_BEAGLE_redshift = 1

# SAME + mass_BEAGLE_stellar
option_BEAGLE_mass = 1

# Changed limit on BEAGLE mass to account for mass offset
option_BEAGLE_mass_lower = 1

# Adjusting upper limit to account for mass offset
option_BEAGLE_mass_upper = 1

# SAME + sfr_BEAGLE_instant
option_BEAGLE_sfr = 1

# all fields
option_all_fields = 0 # otherwise default is clusters

# which redshift set to copy from santini (affects redshift + upper/lower mass cuts + santini MS plot)
#san_bins = [0, 1, 2, 3, 4] # 1.65, 2.5, 3.5, 4.5, 5.5
san_bins = [3]

# minimum chi2 from BEAGLE fit
option_chi2_cut = 1
max_chi2 = 9.5

# cut by median instantaneous BEAGLE sfr 
option_BEAGLE_sfr_cut = 1
min_BEAGLE_sfr_cut = -2.0
max_BEAGLE_sfr_cut = 10.0 # max is 5.72931

# cut by beta - UV slope from Santini
option_beta_cut = 0
min_beta_cut = -3.5
max_beta_cut = 1.0

# PLOTS
plot_BEAGLE_vs_AD = 0
plot_input_to_2sigma = 1
plot_input_to_2sigma_color = 0
plot_input_to_2sigma_heatplot = 0 # HAS TO USE BEAGLE VALUES OF STELLAR MASS AND INSTANT SFR
plot_2sigma_iterations = 0
plot_2sigma_iterations_rf = 0 # rising falling
plot_MS = 1
plot_MS_zoom = 0
plot_bootstrap = 0
print_bootstrap = 0

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

#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_lester_2600_median.npy'
#    sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini_beta_lester_test.npy'    
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
                    'amp_GMM':              np.load(sbf+'amp_GMM.npy')

                    }



    # =============================================================================
    # clusters only?
    # =============================================================================
    if option_all_fields == 1:
        idx1 = (AD['field']%1.0==0.0) # all
    elif option_all_fields == 0:
        idx1 = (AD['field']%2.0==0.0) # clusters
    
    
    # =============================================================================
    # H band mag cut of < 27.5
    # =============================================================================
    idx2 = (AD['H160']<27.5)
    idx = np.logical_and(idx1,idx2)
    
    # =============================================================================
    # Photo-z from median of 6 methods
    # =============================================================================
    zLow = [1.3, 2.0, 3.0, 4.0, 5.0]
    zHigh = [2.0, 3.0, 4.0, 5.0, 6.0]
    zLow = zLow[san_bin]
    zHigh = zHigh[san_bin]
    wLim = (zHigh - zLow) / 2.0
    zLim = zLow + wLim
    idx3 = (abs(AD['ZBEST']-zLim) < wLim)
    idx = np.logical_and(idx,idx3)
    
    # =============================================================================
    # Need sample to be complete above given mass (magnification NOT already included)
    # =============================================================================
    massLow = [8.3, 8.5, 8.8, 8.8, 8.8]
    massLow = massLow[san_bin]
    massLow_offset = np.mean(((data['mass_AD_neb'] + np.log10(AD['MAGNIF'])) - (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF'])))[np.isfinite((data['mass_AD_neb'] + np.log10(AD['MAGNIF'])) - (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF'])))]) # (ignore nan values)
    
    if option_BEAGLE_mass == 0:
        idx4 = (data['mass_AD_neb'] + np.log10(AD['MAGNIF']) > massLow)
        idx = np.logical_and(idx,idx4)
    elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_lower== 0:
        # changing mass to BEAGLE
        idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > massLow)
        idx = np.logical_and(idx,idx4)
    elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_lower == 1:
        # adjusting for offset in mass between AD and BEAGLE (ignore nan values)
        idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > massLow-massLow_offset)
        idx = np.logical_and(idx,idx4)
    
    # =============================================================================
    # Adding RELFLAG just to be sure...
    # =============================================================================
    idx5 = (AD['RELFLAG']==1.0)
    idx = np.logical_and(idx,idx5)
    
    # =============================================================================
    # Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
    # =============================================================================
    massHigh = [10.2, 10.6, 10.8, 10.8, 10.8]
    massHigh = massHigh[san_bin]
    massHigh_offset = np.mean((data['mass_AD_neb'] - data['mass_BEAGLE_stellar'])[np.isfinite(data['mass_AD_neb'] - data['mass_BEAGLE_stellar'])]) # (ignore nan values)
    
    if option_BEAGLE_mass == 0:
        idx6 = (data['mass_AD_neb'] < massHigh)
        idx = np.logical_and(idx,idx6)
    elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_upper == 0:
        # changing mass to BEAGLE
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh)
        idx = np.logical_and(idx,idx6)    
    elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_upper == 1:
        # adjusting for offset in mass between AD and BEAGLE 
        idx6 = (data['mass_BEAGLE_stellar'] < massHigh-massHigh_offset)
        idx = np.logical_and(idx,idx6)
        
    # =============================================================================
    # BEAGLE REDSHIFT AS WELL AS ASTRODEEP
    # =============================================================================
    if option_BEAGLE_redshift == 1:
        idx7 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
        idx = np.logical_and(idx,idx7)
    
    # =============================================================================
    # minimum chi2 cut from BEAGLE fits
    # =============================================================================
    if option_chi2_cut == 1:
        idx8 = (data['min_chi2_BEAGLE']>0) & (data['min_chi2_BEAGLE']<max_chi2)
        idx = np.logical_and(idx,idx8)
    
    # =============================================================================
    # cut by median instantaneous BEAGLE sfr 
    # =============================================================================
    if option_BEAGLE_sfr_cut == 1:
        idx9 = (data['sfr_BEAGLE_instant']>min_BEAGLE_sfr_cut) & (data['sfr_BEAGLE_instant']<max_BEAGLE_sfr_cut)
        idx = np.logical_and(idx,idx9)

    # =============================================================================
    # cut by beta slope when calculating santini SFR
    # =============================================================================
    if option_beta_cut == 1:
        idx10 = (data['sfr_SAN_beta']>min_beta_cut) & (data['sfr_SAN_beta']<max_beta_cut)
        idx = np.logical_and(idx,idx10)

    # =============================================================================
    # COMBINE THE ABOVE
    # =============================================================================
    
    #for key in data.keys():
    #    print(key, len(data[key]), min(data[key]))
    #    print(data[key])
    
    for key in data.keys():
        data[key] = data[key][idx]
    #    print(key, len(data[key]), min(data[key]))

#    pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data.p','w'))
    
    print(len(data['sfr_SAN_beta']))
    print(sum(idx1), sum(idx5), sum(idx2), sum(idx4), sum(idx3), sum(idx6), sum(idx8), sum(idx9))
    
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
    # CHECK THIS
    alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
    beta_err_san = A_err_san
    
    # =============================================================================
    # FEW PLOTS
    # =============================================================================
    x = np.linspace(7, 10.5)
    
    if plot_BEAGLE_vs_AD == 1:
        plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], alpha=0.1)
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlim(1.3, 2)
        plt.ylim(1.3, 2)
        plt.plot((0,10),(0,10),color='k')
        plt.xlabel('redshift AD')
        plt.ylabel('redshift BEAGLE')
        plt.show()
        
        plt.scatter(data['mass_AD_neb'], data['mass_BEAGLE_stellar'], alpha=0.1, c=data['mag_AD'])
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlim(7, 11)
        plt.ylim(7, 11)
        plt.plot((5,12),(5,12),color='k')
        plt.xlabel('mass AD neb')
        plt.ylabel('mass BEAGLE stellar')
        plt.colorbar()
        plt.show()

#        plt.scatter(data['mass_AD'], data['mass_BEAGLE_stellar'], alpha=0.1, c=data['mag_AD'])
#        plt.title('${} < z < {}$'.format(zLow, zHigh))
#        plt.xlim(7, 11)
#        plt.ylim(7, 11)
#        plt.plot((5,12),(5,12),color='k')
#        plt.xlabel('mass AD')
#        plt.ylabel('mass BEAGLE stellar')
#        plt.colorbar()
#        plt.show()
        
        plt.scatter(data['sfr_SAN'], data['sfr_BEAGLE_instant'], alpha=0.1)
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        plt.plot((-3,3),(-3,3),color='k')
        plt.xlabel('sfr SAN')
        plt.ylabel('sfr BEAGLE instant')
        plt.show()
    
    # =============================================================================
    # 2 sigma process
    # =============================================================================
    
    #color = data['msa_BEAGLE'] # colorbar for iterative 2 sigma clipping plots
    #color = data['tau_BEAGLE']
    #color = data['tauv_BEAGLE']
    color = data['msa_BEAGLE'] - data['tau_BEAGLE']
    
    # color for rising falling
    color_rf = data['msa_BEAGLE'] - data['tau_BEAGLE']
      
    if option_BEAGLE_mass == 0 and option_BEAGLE_sfr == 0:
        idx_data = (data['sfr_SAN']>-90.0) # -99.0 is assigned when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)
        mass = data['mass_AD_neb'][data['sfr_SAN']>-90.0] 
        sfr = data['sfr_SAN'][data['sfr_SAN']>-90.0]
        color = color[data['sfr_SAN']>-90.0]
        color_rf = color_rf[data['sfr_SAN']>-90.0]
    elif option_BEAGLE_mass == 1 and option_BEAGLE_sfr == 0:
        # SAME + mass_BEAGLE_stellar
        idx_data = (data['sfr_SAN']>-90.0)
        mass = data['mass_BEAGLE_stellar'][data['sfr_SAN']>-90.0]
        sfr = data['sfr_SAN'][data['sfr_SAN']>-90.0]
        color = color[data['sfr_SAN']>-90.0]
        color_rf = color_rf[data['sfr_SAN']>-90.0]
    elif option_BEAGLE_mass == 0 and option_BEAGLE_sfr == 1:
        idx_data = (data['sfr_BEAGLE_instant']>-30.0)
        mass = data['mass_AD_neb'][data['sfr_BEAGLE_instant']>-30.0] # -101 AD OBJECT WAS NOT A BEAGLE INPUT, -102 AD OBJECT WAS NOT FITTED BY BEAGLE, -30 was during creation of instant sfr
        sfr = data['sfr_BEAGLE_instant'][data['sfr_BEAGLE_instant']>-30.0]
        color = color[data['sfr_BEAGLE_instant']>-30.0]
        color_rf = color_rf[data['sfr_BEAGLE_instant']>-30.0]
    elif option_BEAGLE_mass == 1 and option_BEAGLE_sfr == 1:   
        # SAME + sfr_BEAGLE_instant
        idx_data = (data['sfr_BEAGLE_instant']>-30.0)
        mass = data['mass_BEAGLE_stellar'][data['sfr_BEAGLE_instant']>-30.0]
        sfr = data['sfr_BEAGLE_instant'][data['sfr_BEAGLE_instant']>-30.0]
        color = color[data['sfr_BEAGLE_instant']>-30.0]
        color_rf = color_rf[data['sfr_BEAGLE_instant']>-30.0]  
    
    for key in data.keys():
        data[key] = data[key][idx_data]   
    
    outliers = 1
    
    # values saved as _0 prior to 2sigma clipping
    mass0 = mass
    sfr0 = sfr
    fit0 = np.polyfit(mass0, sfr0, 1)
    sfr_residuals0 = sfr0 - (fit0[0]*mass0 + fit0[1])
    sigma0 = np.std(sfr_residuals0)
    
    print(len(data['sfr_SAN_beta']))
    
    # =============================================================================
    # some plots
    # =============================================================================
    if plot_input_to_2sigma == 1:
        plt.scatter(mass0, sfr0)
        plt.show()
    
    if plot_input_to_2sigma_color == 1:
        plt.figure(figsize=(10, 6))    
        plt.scatter(mass0, sfr0, alpha=0.3, c=color)
        plt.xlim(7, 10)
        plt.ylim(-4, 3)
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.colorbar()
        plt.show()
        print('BEFORE 2sig LEN MASS', len(mass))
    
    
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
        color = color[idx]
        color_rf = color_rf[idx]
        
        # =============================================================================
        # PLOT ITERATIONS - colour coded
        # =============================================================================
        if plot_2sigma_iterations == 1:
            plt.scatter(mass, sfr, alpha=0.3, c=color)
            plt.xlim(7, 10)
            plt.ylim(-4, 3)
            plt.colorbar()
            plt.show()
            print(fit, sigma)   
            print('LEN MASS', len(mass))
    
        # =============================================================================
        # PLOT ITERATIONS - rising falling - color must be correctly set
        # =============================================================================
        if plot_2sigma_iterations_rf == 1:
            idx_falling = (color_rf > 0)
            idx_rising = (color_rf < 0)
            plt.scatter(mass[idx_rising], sfr[idx_rising], alpha=0.3, color='k')
            plt.scatter(mass[idx_falling], sfr[idx_falling], alpha=0.3, color='r')
            plt.xlim(7, 10)
            plt.ylim(-4, 3)
            plt.show()
            print(fit, sigma)   
            print('LEN MASS', len(mass))
        
    fit_clipping = fit

    
    # =============================================================================
    # FIRST YEAR REPORT
    # =============================================================================
    if plot_MS == 1:
        plt.figure(figsize=(6, 6))
        plt.title('${} < z < {}$'.format(zLow, zHigh))
        plt.xlabel(r'$M_\mathrm{tot}$')
        plt.ylabel(r'$\Psi$')
        plt.scatter(mass0, sfr0, marker='x')
        plt.scatter(mass, sfr, marker='x')
        plt.plot(x, fit0[0]*x + fit0[1], color='#1f77b4', label='Without clipping', linewidth=2)
        plt.plot(x, fit[0]*x + fit[1], color='#d62728', label='With clipping', linewidth=2)
        plt.plot(x, beta_san[san_bin]*x + alpha_san[san_bin], color='#2ca02c', label='Santini+17', linewidth=2)
                 
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
    
    iterations = 10 # 10000 for report
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




