#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 10:37:20 2021

@author: lester
"""

# https://towardsdatascience.com/an-introduction-to-making-scientific-publication-plots-with-python-ea19dfa7f51e

# Import required packages
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import matplotlib.font_manager as fm
import pickle
from scipy.stats import norm, multivariate_normal
import corner
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
import copy
import sys

# Collect all the font names available to matplotlib
#font_names = [f.name for f in fm.fontManager.ttflist]
#print(font_names)

# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
mpl.rc('image', cmap='jet')
cmap = mpl.cm.get_cmap('jet')
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['text.usetex'] = True
fontsize_legend = 10
fontsize_axes = 24
figuresize = 7

normalisation = 9.7

load = True
# load = False
save = False
 
# =============================================================================
# useful axis strings
# =============================================================================

string_slope = r'$\mathrm{MS\,Slope,}\,\beta$'
string_normalisation = r'MS Normalisation, $\alpha_{9.7}$'
string_scatter = r'$\mathrm{MS\,Scatter,}\,\sigma$'
string_ssfr = r'$\log(\mathrm{sSFR} \, / \, \mathrm{yr}^{-1})$'
string_mass = r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$'
string_sfr = r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$'
string_deltaMS = r'$\Delta_{MS}$'
string_prob_ratio = r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$'
string_bias_test = r'$\Delta \mathrm{Parameter}$'
string_pbad = r'pbad'
string_outlier_mean = r'outlier mean'
string_outlier_sigma = r'outlier sigma'



# =============================================================================
# faff
# =============================================================================
'''
    
    ssfr_a_arr = np.array([chain_MS_29_c1['alphaN_a']]).T
    ssfr_b_arr = np.array([chain_MS_29_c1['alphaN_b']]).T
    ssfr_arr = np.log10(ssfr_a_arr*(1.0+0.0)**ssfr_b_arr) - 9.0
    alpha_arr = ssfr_arr + normalisation
    print(ssfr_a_arr)
    
    plt.figure(figsize=(10, 10))
    plt.plot(alpha_arr)
    plt.show()
    
    ssfr_a_arr = np.array([chain_MS_29_c1k['alphaN_a']]).T
    ssfr_b_arr = np.array([chain_MS_29_c1k['alphaN_b']]).T
    ssfr_arr = np.log10(ssfr_a_arr*(1.0+0.0)**ssfr_b_arr) - 9.0
    alpha_arr = ssfr_arr + normalisation
    #print(ssfr_a_arr)
    
    plt.figure(figsize=(10, 10))
    plt.plot(alpha_arr)
    plt.show()
'''



# =============================================================================
# opening data
# =============================================================================
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3

def read_chain(zlow, zhigh, chain_MS, fit):
    
    z = np.linspace(zlow, zhigh, 1000)

    redshift_arr = np.repeat(np.array([z]).T, len(chain_MS), axis=1).T

    beta_a_arr = np.array([chain_MS['beta_a']]).T
    beta_b_arr = np.array([chain_MS['beta_b']]).T
    beta_arr = beta_a_arr + redshift_arr*beta_b_arr
    beta_16_arr = np.percentile(beta_arr, 16, axis=0)
    beta_50_arr = np.median(beta_arr, axis=0)
    beta_84_arr = np.percentile(beta_arr, 84, axis=0)

    if fit == 'const_beta_ssfr_alpha' or fit == 'const_beta_ssfr_alpha_fixed_pbad':

        ssfr_a_arr = np.array([chain_MS['alphaN_a']]).T
        ssfr_b_arr = np.array([chain_MS['alphaN_b']]).T
        ssfr_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr) - 9.0
        ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
        ssfr_50_arr = np.median(ssfr_arr, axis=0)
        ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)
        
        alpha_arr = ssfr_arr + normalisation
        alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
        alpha_50_arr = np.median(alpha_arr, axis=0)
        alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)
        
    elif fit == 'const_beta_linear_alpha':

        alpha_a_arr = np.array([chain_MS['alphaN_a']]).T
        alpha_b_arr = np.array([chain_MS['alphaN_b']]).T        
        alpha_arr = alpha_a_arr + redshift_arr*alpha_b_arr
        alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
        alpha_50_arr = np.median(alpha_arr, axis=0)
        alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)

        ssfr_arr = alpha_arr - normalisation
        ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
        ssfr_50_arr = np.median(ssfr_arr, axis=0)
        ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)
        
    if fit == 'const_beta_ssfr_alpha_fixed_pbad':
        pbad_arr = np.array([chain_MS['pbad'][:,0]]).T + redshift_arr*0
        pbad_16_arr = np.percentile(pbad_arr, 16, axis=0)
        pbad_50_arr = np.median(pbad_arr, axis=0)
        pbad_84_arr = np.percentile(pbad_arr, 84, axis=0)
    else:
        pbad_arr = np.array([chain_MS['pbad']]).T + redshift_arr*0
        pbad_16_arr = np.percentile(pbad_arr, 16, axis=0)
        pbad_50_arr = np.median(pbad_arr, axis=0)
        pbad_84_arr = np.percentile(pbad_arr, 84, axis=0)
          
#        print(chain_MS['pbad'][:,0])
#        chain_MS['pbad'] = chain_MS['pbad'][:,0]
#        exit
        
    sig0_arr = np.array([chain_MS['sig0']]).T + redshift_arr*0
    sig0_16_arr = np.percentile(sig0_arr, 16, axis=0)
    sig0_50_arr = np.median(sig0_arr, axis=0)
    sig0_84_arr = np.percentile(sig0_arr, 84, axis=0)

    k_arr = np.array([chain_MS['k']]).T + redshift_arr*0
    k_16_arr = np.percentile(k_arr, 16, axis=0)
    k_50_arr = np.median(k_arr, axis=0)
    k_84_arr = np.percentile(k_arr, 84, axis=0)
    
    outlier_mean_arr = np.array([chain_MS['outlier_mean']]).T + redshift_arr*0
    outlier_mean_16_arr = np.percentile(outlier_mean_arr, 16, axis=0)
    outlier_mean_50_arr = np.median(outlier_mean_arr, axis=0)
    outlier_mean_84_arr = np.percentile(outlier_mean_arr, 84, axis=0)
    
    outlier_sigma_arr = np.array([chain_MS['outlier_sigma']]).T + redshift_arr*0
    outlier_sigma_16_arr = np.percentile(outlier_sigma_arr, 16, axis=0)
    outlier_sigma_50_arr = np.median(outlier_sigma_arr, axis=0)
    outlier_sigma_84_arr = np.percentile(outlier_sigma_arr, 84, axis=0)
    
    # DICTIONARY
    
    dic = {}
    
    dic['z'] = z # array of 1000 numbers between zlow and zhigh
    
    dic['beta_16_arr'] = beta_16_arr # 1000 values (1 per redshift interval)
    dic['beta_50_arr'] = beta_50_arr # 1000 values (1 per redshift interval)
    dic['beta_84_arr'] = beta_84_arr # 1000 values (1 per redshift interval)
    
    dic['ssfr_16_arr'] = ssfr_16_arr # 1000 values (1 per redshift interval)
    dic['ssfr_50_arr'] = ssfr_50_arr # 1000 values (1 per redshift interval)
    dic['ssfr_84_arr'] = ssfr_84_arr # 1000 values (1 per redshift interval)

    dic['alpha_16_arr'] = alpha_16_arr # 1000 values (1 per redshift interval)
    dic['alpha_50_arr'] = alpha_50_arr # 1000 values (1 per redshift interval)
    dic['alpha_84_arr'] = alpha_84_arr # 1000 values (1 per redshift interval)
    
    dic['sig0_16_arr'] = sig0_16_arr # 1000 values (1 per redshift interval)
    dic['sig0_50_arr'] = sig0_50_arr # 1000 values (1 per redshift interval)
    dic['sig0_84_arr'] = sig0_84_arr # 1000 values (1 per redshift interval)

    dic['k_16_arr'] = k_16_arr # 1000 values (1 per redshift interval)
    dic['k_50_arr'] = k_50_arr # 1000 values (1 per redshift interval)
    dic['k_84_arr'] = k_84_arr # 1000 values (1 per redshift interval)
    
    dic['pbad_16_arr'] = pbad_16_arr # 1000 values (1 per redshift interval)
    dic['pbad_50_arr'] = pbad_50_arr # 1000 values (1 per redshift interval)
    dic['pbad_84_arr'] = pbad_84_arr # 1000 values (1 per redshift interval)

    dic['outlier_mean_16_arr'] = outlier_mean_16_arr # 1000 values (1 per redshift interval)
    dic['outlier_mean_50_arr'] = outlier_mean_50_arr # 1000 values (1 per redshift interval)
    dic['outlier_mean_84_arr'] = outlier_mean_84_arr # 1000 values (1 per redshift interval)    
    
    dic['outlier_sigma_16_arr'] = outlier_sigma_16_arr # 1000 values (1 per redshift interval)
    dic['outlier_sigma_50_arr'] = outlier_sigma_50_arr # 1000 values (1 per redshift interval)
    dic['outlier_sigma_84_arr'] = outlier_sigma_84_arr # 1000 values (1 per redshift interval)    
    
    return dic

#%%
# =============================================================================
# get medians for MS plots
# =============================================================================

def get_medians(chain_MS):
#    names = chain_original.dtype.names
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
    dic = {}
    for name in names:
        dic[name] = np.median(chain_MS[name])
    return dic

if load:

#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
#        chain_MS_29_c1 = pickle.load(f, encoding='latin1')  
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z2p0-3p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
#        chain_MS_29_c2 = pickle.load(f, encoding='latin1')  
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z3p0-4p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
#        chain_MS_29_c3 = pickle.load(f, encoding='latin1')  
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z4p0-5p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
#        chain_MS_29_c4 = pickle.load(f, encoding='latin1')  
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z5p0-6p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
#        chain_MS_29_c5 = pickle.load(f, encoding='latin1') 
#    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha.p', 'rb') as f:
#        chain_MS_29_c = pickle.load(f, encoding='latin1')
#    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_k_fitted_1_2.p', 'rb') as f:
#        chain_MS_29_c1k = pickle.load(f, encoding='latin1')
#        


    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000.p', 'rb') as f:
        chain_MS_29_c1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z2p0-3p0_4x5000.p', 'rb') as f:
        chain_MS_29_c2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z3p0-4p0_4x5000.p', 'rb') as f:
        chain_MS_29_c3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z4p0-5p0_4x5000.p', 'rb') as f:
        chain_MS_29_c4 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z5p0-6p0_4x5000.p', 'rb') as f:
        chain_MS_29_c5 = pickle.load(f, encoding='latin1') 
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha.p', 'rb') as f:
        chain_MS_29_c = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_k_fitted_1_2.p', 'rb') as f:
        chain_MS_29_c1k = pickle.load(f, encoding='latin1')
        

#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z1p25-2p0_4x5000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_mock.p', 'rb') as f:
#        chain_MS_mock_z1 = pickle.load(f, encoding='latin1')          
    '''    
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm8p5.p', 'rb') as f:
        chain_MS_29_c1_lm8p5 = pickle.load(f, encoding='latin1')  
        
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm9p0.p', 'rb') as f:
        chain_MS_29_c1_lm9p0 = pickle.load(f, encoding='latin1')  
        '''
    # AD catalogue
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)

    # from /Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/replicating_santini_with_santini_input_3dGMM.py
    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p', 'rb') as f:
        ADx = pickle.load(f, encoding='latin1') 
#    print(ADx.keys())        

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p90_new.fits'
    mass_completeness_limits_0p90_new = fits.open(fileName)
    #print(data_fits.info())
    #print(data_fits[1].header)

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z1p25-2p0.fits'
    s29z1 = fits.open(fileName)[1].data 
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z2p0-3p0.fits'
    s29z2 = fits.open(fileName)[1].data 
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z3p0-4p0.fits'
    s29z3 = fits.open(fileName)[1].data 
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z4p0-5p0.fits'
    s29z4 = fits.open(fileName)[1].data 
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z5p0-6p0.fits'
    s29z5 = fits.open(fileName)[1].data 
    
#    rc_c1 = read_chain(1.25, 2.0, chain_MS_29_c1, 'const_beta_linear_alpha')
#    rc_c1_lm8p5 = read_chain(1.25, 2.0, chain_MS_29_c1_lm8p5, 'const_beta_linear_alpha')
#    rc_c1_lm9p0 = read_chain(1.25, 2.0, chain_MS_29_c1_lm9p0, 'const_beta_linear_alpha')
#    rc_c2 = read_chain(2.0, 3.0, chain_MS_29_c2, 'const_beta_linear_alpha')
#    rc_c3 = read_chain(3.0, 4.0, chain_MS_29_c3, 'const_beta_linear_alpha')
#    rc_c4 = read_chain(4.0, 5.0, chain_MS_29_c4, 'const_beta_linear_alpha')
#    rc_c5 = read_chain(5.0, 6.0, chain_MS_29_c5, 'const_beta_linear_alpha')
#    rc_c = read_chain(1.25, 6.0, chain_MS_29_c, 'const_beta_ssfr_alpha')
#    rc_c1k = read_chain(1.25, 2.0, chain_MS_29_c1k, 'const_beta_ssfr_alpha')
#    rc_mock_z1 = read_chain(1.25, 2.0, chain_MS_mock_z1, 'const_beta_linear_alpha')
    

    # FOR POSTER (from original red line in KICC report)
    rc_c1 = read_chain(1.25, 2.0, chain_MS_29_c1, 'const_beta_ssfr_alpha')
    rc_c2 = read_chain(2.0, 3.0, chain_MS_29_c2, 'const_beta_ssfr_alpha')
    rc_c3 = read_chain(3.0, 4.0, chain_MS_29_c3, 'const_beta_ssfr_alpha')
    rc_c4 = read_chain(4.0, 5.0, chain_MS_29_c4, 'const_beta_ssfr_alpha')
    rc_c5 = read_chain(5.0, 6.0, chain_MS_29_c5, 'const_beta_ssfr_alpha')
    rc_c = read_chain(1.25, 6.0, chain_MS_29_c, 'const_beta_ssfr_alpha')
    rc_c1k = read_chain(1.25, 2.0, chain_MS_29_c1k, 'const_beta_ssfr_alpha')

#%%

#%%
# =============================================================================
# redshift bin MS colour coded by MS probability SINGLE PLOT
# =============================================================================

s = [s29z1, s29z2, s29z3, s29z4, s29z5]
rc = [rc_c1, rc_c2, rc_c3, rc_c4, rc_c5]
z_lower = [1.25, 2.0, 3.0, 4.0, 5.0]
z_upper = [2.0, 3.0, 4.0, 5.0, 6.0]
medians = [get_medians(chain_MS_29_c1), get_medians(chain_MS_29_c2), get_medians(chain_MS_29_c3), get_medians(chain_MS_29_c4), get_medians(chain_MS_29_c5)]


fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
ax2 = fig.add_axes([0.85, 0, 0.85, 0.84]) #[left, bottom, width, height]
ax3 = fig.add_axes([1.70, 0, 0.85, 0.84]) #[left, bottom, width, height]
ax4 = fig.add_axes([2.55, 0, 0.85, 0.84]) #[left, bottom, width, height]
ax5 = fig.add_axes([3.40, 0, 0.85, 0.84]) #[left, bottom, width, height]

ax = [ax1, ax2, ax3, ax4, ax5]

for i in range(len(s)):

    ### ssfr 
    sfr_surface_real = ((medians[i]['beta_a']+s[i]['redshift_BEAGLE']*medians[i]['beta_b'])*s[i]['mass_BEAGLE_stellar'])+(np.log10(medians[i]['alphaN_a']*((1+s[i]['redshift_BEAGLE'])**medians[i]['alphaN_b'])))+normalisation-9.0-(normalisation*(medians[i]['beta_a']+s[i]['redshift_BEAGLE']*medians[i]['beta_b']))
    
    ### linear
    #sfr_surface_real = ((medians['beta_a']+s['redshift_BEAGLE']*medians['beta_b'])*(s['mass_BEAGLE_stellar'] - 9.7))  + \
    #                    (medians['alphaN_a']+s['redshift_BEAGLE']*medians['alphaN_b'])
    
    log_p_xi_eta_theta = norm.logpdf(s[i]['sfr_BEAGLE_instant'], scale=medians[i]['sig0'], loc=sfr_surface_real)
    log_p_eta_xi_theta = norm.logpdf(s[i]['sfr_BEAGLE_instant'], scale=medians[i]['sig0'], loc=sfr_surface_real)
    p_bad = norm.pdf(s[i]['sfr_BEAGLE_instant'], scale=medians[i]['outlier_sigma'], loc=medians[i]['outlier_mean'])
    
    z_bad = medians[i]['pbad']*p_bad
    z_good = (1.0-medians[i]['pbad'])*np.exp(log_p_eta_xi_theta)
    
    idx_sort = np.argsort(z_good/z_bad)
    
    ax[i].plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')
    
    scatter = ax[i].scatter(s[i]['mass_BEAGLE_stellar'][idx_sort], s[i]['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5, s=100)
    
    print(len(s[i]['mass_BEAGLE_stellar']))
    
    xlow = 6.5
    xhigh = 10.5
    ylow = -2.5
    yhigh = 2.5
    
    xi_min = 8.5
    xi_max = 10.0
    
    lw = 3
    
    ax[i].plot((9.7,9.7), (ylow, yhigh), color='gray', linestyle='dashed', linewidth=2)
    
    x_tmp = np.array([xlow, xhigh])
    ax[i].plot(x_tmp, (x_tmp-9.7)*rc[i]['beta_50_arr'][0] + rc[i]['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, {} $<$ z $<$ {}'.format(str(z_lower[i]), str(z_upper[i])))
     

    sig = rc[i]['sig0_50_arr'][0] * ( ((1.0-rc[i]['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
    ax[i].plot(x_tmp, (x_tmp-9.7)*rc[i]['beta_50_arr'][0] + rc[i]['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
    ax[i].plot(x_tmp, (x_tmp-9.7)*rc[i]['beta_50_arr'][0] + rc[i]['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)
    
    # SANTINI, log(SFR) = alpha log(M / M_9p7) + beta
    # 1.3 < z < 2
    s1_alpha = 1.04
    s1_beta = 1.01
    s1_alpha_err = 0.03
    s1_beta_err = 0.04
    
    # santini z=1
    s_x = np.array([8.4, 9.2, 10.0])
    s_y = s1_alpha*(s_x - 9.7) + s1_beta
    s_y_intrinsic = np.array([0.36, 0.35, 0.0])
    s_y_observed = np.array([0.51, 0.46, 0.26])
    
    #ax1.plot(s_x, s_y, color='k', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    #ax1.plot(s_x, s_y+s_y_intrinsic, color='k', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    #ax1.plot(s_x, s_y-s_y_intrinsic, color='k', linestyle='--', linewidth=lw)
    #ax1.plot(s_x, s_y+s_y_observed, color='k', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    #ax1.plot(s_x, s_y-s_y_observed, color='k', linestyle=':', linewidth=lw)
    
    ax[i].set_xlim(xlow, xhigh)
    ax[i].set_xlabel(r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
    ax[i].xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax[i].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax[i].xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax[i].xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax[i].xaxis.set_tick_params(labelsize=fontsize_axes)
    ax[i].legend(loc='lower right', frameon=True, fontsize=1.44*fontsize_legend, framealpha=1.0)
    
    ax[i].text(6.7, 2.0, r'{} $<$ z $<$ {}'.format(str(z_lower[i]), str(z_upper[i])))
    ax[i].set_ylim(ylow, yhigh)
    
    if i != 0:
        ax[i].yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
        ax[i].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
        ax[i].yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
        ax[i].yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
        ax[i].yaxis.set_tick_params(labelsize=0*fontsize_axes)
    
    else:
        ax[i].set_ylabel(r'$\log(\mathrm{SFR} \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
        ax[i].yaxis.set_tick_params(labelsize=fontsize_axes)

cbaxes = fig.add_axes([4.25, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\mathrm{MS}\longleftarrow \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \longrightarrow \mathrm{Outliers}$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=0, width=2, direction='in', labelsize=0*fontsize_axes)

if save:
    plt.savefig('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/EAS_poster/1.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()



#%%
# =============================================================================
# ssfr and scatter vs redshift
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1.5*figuresize))
xlow = -0.3
xhigh = 7.3

param = ['ssfr']

ax1 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1]

for ax in [ax1]:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.xaxis.set_tick_params(which='minor', size=0, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=0, width=2, direction='out')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

for rc in [rc_c]:
    for a, ax in enumerate([ax1]):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='red', linewidth=lw)
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='red', zorder=0)
     
# layout
ylim_low = [-10.5]
ylim_high = [-7.3]
string = [string_ssfr]
ticker_maj = [1.0]
ticker_min = [0.1, 0.1]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

# literature ssfr
mpl.rcParams.update({'errorbar.capsize': 5})

# Salim+07
z_s07 = np.array([0.0992])
ssfr_s07 = np.array([-0.8106]) - 9.0
ssfr_p_s07 = np.array([-0.296]) - ssfr_s07 - 9.0
ssfr_m_s07 = ssfr_s07 + 9.0 - np.array([-1.2918])

# Santini+09
z_s09 = np.array([0.4429, 0.8013, 1.2547, 2.0007])
ssfr_s09 = np.array([-0.2388, -0.0578, 0.1185, 0.2995]) - 9.0
ssfr_p_s09 = np.array([0.0566, 0.2376, 0.4187, 0.595]) - ssfr_s09 - 9.0
ssfr_m_s09 = ssfr_s09 + 9.0 - np.array([-0.5342, -0.358, -0.1817, 4.1691e-3])

# Santini+17 Observed
z_s17 = np.array([1.657, 2.5127, 3.5072, 4.502, 5.5112])
ssfr_s17 = np.array([0.2662, 0.4282, 0.4949, 0.2424, 0.5615]) - 9.0
ssfr_p_s17 = np.array([0.5855, 0.5855, 0.7379, 0.4235, 1.2144]) - ssfr_s17 - 9.0
ssfr_m_s17 = ssfr_s17 + 9.0 - np.array([-0.0482, 0.2852, 0.2472, 0.0661, -0.0911])

# Santini+17 Simulation
z_s17s = np.array([1.657, 2.498, 3.5, 4.5093, 5.5039])
ssfr_s17s = np.array([0.3139, 0.5188, 0.6712, 0.6665, 1.2953]) - 9.0

# Reddy+12
z_r12 = np.array([2.3006, 3.0027])
ssfr_r12 = np.array([0.3663, 0.3472]) - 9.0
ssfr_p_r12 = np.array([0.7475, 0.7236]) - ssfr_r12 - 9.0
ssfr_m_r12 = ssfr_r12 + 9.0 - np.array([-5.3603e-3, -0.0244])

# de Barros+14
z_b14 = np.array([3.2294, 3.7121, 4.8237, 5.9062])
ssfr_b14 = np.array([0.8046, 0.8046, 1.2382, 1.2335]) - 9.0
ssfr_p_b14 = np.array([1.3002, 1.3145, 1.3096, 1.2906]) - ssfr_b14 - 9.0
ssfr_m_b14 = ssfr_b14 + 9.0 - np.array([-0.0578, -0.0482, 0.0899, 0.0995])

# Stark+13
z_s13 = np.array([3.7999, 4.9774, 5.8916, 5.8916, 6.7911, 6.7911])
ssfr_s13 = np.array([0.7522, 0.7332, 0.8332, 0.7808, 1.1191, 0.9524]) - 9.0
ssfr_p_s13 = np.array([1.0715, 1.0476, 1.1048, 1.1048, 1.2859, 1.2859]) - ssfr_s13 - 9.0
ssfr_m_s13 = ssfr_s13 + 9.0 - np.array([0.4663, 0.4378, 0.4997, 0.4997, 0.676, 0.676])

# Gonzalez+14
z_g14 = np.array([3.8072, 5.0212, 5.9135])
ssfr_g14 = np.array([0.5426, 0.5331, 0.6808]) - 9.0
ssfr_p_g14 = np.array([0.8285, 0.8237, 0.9762]) - ssfr_g14 - 9.0
ssfr_m_g14 = ssfr_g14 + 9.0 - np.array([0.2376, 0.2328, 0.3806])

# Bouwens+12
z_b12 = np.array([4.0046, 5.0066, 6.0086, 7.208])
ssfr_b12 = np.array([0.7141, 0.6855, 0.4997, 0.7284]) - 9.0
ssfr_p_b12 = np.array([1.0048, 0.9762, 0.7904, 0.8475]) - ssfr_b12 - 9.0
ssfr_m_b12 = ssfr_b12 + 9.0 - np.array([0.4139, 0.3853, 0.1995, 0.6236])

# Marmol-Queralto+16
z_m16 = np.array([4.385])
ssfr_m16 = np.array([0.7284]) - 9.0
ssfr_p_m16 = np.array([1.0286]) - ssfr_m16 - 9.0
ssfr_m_m16 = ssfr_m16 + 9.0 - np.array([0.433])

# Menci et al. (2014) SAM
z_m14 = np.array([0.1504, 0.3624, 0.5819, 0.8525, 1.0719, 1.4522, 1.8763, 2.403, 2.9149, 3.4487, 4.0632, 4.8969, 5.6429, 6.228, 6.7473])
ssfr_m14 = np.array([-1.0965, -0.9297, -0.8154, -0.6915, -0.5867, -0.458, -0.3389, -0.215, -0.1197, -0.0339, 0.0661, 0.1948, 0.3139, 0.4139, 0.514]) - 9.0


ax1.scatter(z_s07, ssfr_s07, label='Salim+07')
ax1.errorbar(z_s07, ssfr_s07, yerr=(ssfr_m_s07, ssfr_p_s07), ls='none')

ax1.scatter(z_s09, ssfr_s09, label='Santini+09')
ax1.errorbar(z_s09, ssfr_s09, yerr=(ssfr_m_s09, ssfr_p_s09), ls='none')

ax1.scatter(z_s17, ssfr_s17, label='Santini+17 Observed')
ax1.errorbar(z_s17, ssfr_s17, yerr=(ssfr_m_s17, ssfr_p_s17), ls='none')

ax1.scatter(z_s17s, ssfr_s17s, label='Santini+17 Simulation')
ax1.errorbar(0, 0, ls='none')

ax1.scatter(z_r12, ssfr_r12, label='Reddy+12')
ax1.errorbar(z_r12, ssfr_r12, yerr=(ssfr_m_r12, ssfr_p_r12), ls='none')

ax1.scatter(z_b14, ssfr_b14, label='de Barros+14')
ax1.errorbar(z_b14, ssfr_b14, yerr=(ssfr_m_b14, ssfr_p_b14), ls='none')

ax1.scatter(z_s13, ssfr_s13, label='Stark+13')
ax1.errorbar(z_s13, ssfr_s13, yerr=(ssfr_m_s13, ssfr_p_s13), ls='none')

ax1.scatter(z_g14, ssfr_g14, label='Gonzalez+14')
ax1.errorbar(z_g14, ssfr_g14, yerr=(ssfr_m_g14, ssfr_p_g14), ls='none')

ax1.scatter(z_b12, ssfr_b12, label='Bouwens+12')
ax1.errorbar(z_b12, ssfr_b12, yerr=(ssfr_m_b12, ssfr_p_b12), ls='none')

ax1.scatter(z_m16, ssfr_m16, label='Marmol-Queralto+16')
ax1.errorbar(z_m16, ssfr_m16, yerr=(ssfr_m_m16, ssfr_p_m16), ls='none')


ax1.plot(0, 0, color='red', label='Our Work', linewidth=lw)
ax1.plot(z_m14, ssfr_m14, label='Menci et al. (2014) SAM', color='gray', linestyle=':')

z_225 = np.linspace(0, 10, 1000)
idx = np.absolute(rc_c['z']-2.0).argmin()
ax1.plot(z_225, (np.log10((1+z_225)**2.25) - 9.0) * (rc_c['ssfr_50_arr'][idx] / (np.log10((1+2.0)**2.25) - 9.0)), color='gray', linestyle='--', label='$\sim(1+z)^{2.25}$') # normalise Dekel et al. 2009 to z=2

ax1.legend(loc='lower right', frameon=False, fontsize=1.26*fontsize_legend, ncol=2)

plt.text(0, -7.6, r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) \sim 9.7$')

ax1.set_xlabel('Redshift', labelpad=10)
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

if save:
    plt.savefig('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/EAS_poster/2.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# combining 2 and 4
# =============================================================================

fig = plt.figure(figsize=(0.8*figuresize, 0.8*figuresize))
medians = get_medians(chain_MS_29_c1k)

sfr_surface_real = ((medians['beta_a']+s29z1['redshift_BEAGLE']*medians['beta_b'])*s29z1['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+s29z1['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+s29z1['redshift_BEAGLE']*medians['beta_b']))

log_p_xi_eta_theta = norm.logpdf(s29z1['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(s29z1['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
p_bad = norm.pdf(s29z1['sfr_BEAGLE_instant'], scale=medians['outlier_sigma'], loc=medians['outlier_mean'])

z_bad = medians['pbad']*p_bad
z_good = (1-medians['pbad'])*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]


ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

scatter = ax1.scatter(s29z1['mass_BEAGLE_stellar'][idx_sort], s29z1['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5, s=100)

print(len(s29z1['mass_BEAGLE_stellar']))

xlow = 6.5
xhigh = 10.5
ylow = -2.5
yhigh = 2.5

ximin = 8.5
ximax = 10.0

lw = 3

ax1.plot((9.7,9.7), (ylow, yhigh), color='gray', linestyle='dashed', linewidth=2)

x_tmp = np.array([xlow, xhigh])
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, 1.25 $<$ z $<$ 2.0')
 
sig = rc_c1k['sig0_50_arr'][0] * ( ((1.0-rc_c1k['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)

# SANTINI, log(SFR) = alpha log(M / M_9p7) + beta
# 1.3 < z < 2
s1_alpha = 1.04
s1_beta = 1.01
s1_alpha_err = 0.03
s1_beta_err = 0.04

# santini z=1
s_x = np.array([8.4, 9.2, 10.0])
s_y = s1_alpha*(s_x - 9.7) + s1_beta
s_y_intrinsic = np.array([0.36, 0.35, 0.0])
s_y_observed = np.array([0.51, 0.46, 0.26])

ax1.plot(s_x, s_y, color='k', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
ax1.plot(s_x, s_y+s_y_intrinsic, color='k', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
ax1.plot(s_x, s_y-s_y_intrinsic, color='k', linestyle='--', linewidth=lw)
ax1.plot(s_x, s_y+s_y_observed, color='k', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
ax1.plot(s_x, s_y-s_y_observed, color='k', linestyle=':', linewidth=lw)

ax1.set_xlim(xlow, xhigh)
ax1.set_xlabel(r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(ylow, yhigh)
ax1.set_ylabel(r'$\log(\mathrm{SFR} \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

ax1.legend(loc='lower right', frameon=True, fontsize=1.44*fontsize_legend, framealpha=1.0)


cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\mathrm{MS}\longleftarrow \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \longrightarrow \mathrm{Outliers}$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=0, width=2, direction='in', labelsize=0*fontsize_axes)





ax2 = fig.add_axes([1.0, 0, 0.85, 0.84]) #[left, bottom, width, height]

ax2.hist(chain_MS_29_c1k['sig0']*chain_MS_29_c1k['k'], alpha=0.3, bins=30, density=True, range=[0.1, 0.5], label=r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 8.5$')
ax2.hist(chain_MS_29_c1k['sig0'], alpha=0.3, bins=30, density=True, range=[0.1, 0.5], label=r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$')

'''
$\mathrm{Stellar} \,  \mathrm{Mass} = 10^{8.5}\, \mathrm{M_{\odot}}$
$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 8.5$

$\mathrm{Stellar} \,  \mathrm{Mass} = 10^{10.0}\, \mathrm{M_{\odot}}$
$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$
'''

#ax2.hist(((ADx_subset['mass_BEAGLE_stellar']-ADx_subset['mass_SANTINI'])/(1.0+ADx_subset['mass_SANTINI'])), alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
#ax2.hist(((ADx_subset['sfr_BEAGLE_instant']-ADx_subset['sfr_SANTINI'])/(1.0+ADx_subset['sfr_SANTINI'])), alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])

#ax2.hist(ADx_subset['redshift_AD'], alpha=0.3, bins=30)
#ax2.hist(ADx_subset['redshift_BEAGLE'], alpha=0.3, bins=30)

#ax2.set_xlim(-50, 50)
ax2.set_xlabel(r'$\mathrm{Intrinsic} \,  \mathrm{Scatter}$', labelpad=10)
#ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax2.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax2.set_ylim(0.0, 10.0)
#ax2.set_ylabel(r'Count', labelpad=10)
#ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax2.yaxis.set_tick_params(which='major', size=0, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(which='minor', size=0, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(labelsize=0*fontsize_axes)


ax2.legend(loc='lower right', frameon=True, fontsize=1.44*fontsize_legend)




#save = True
if save:
    plt.savefig('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/EAS_poster/3.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()




#%%
# =============================================================================
# just left side of above
# =============================================================================

fig = plt.figure(figsize=(0.8*figuresize, 0.8*figuresize))
medians = get_medians(chain_MS_29_c1k)

sfr_surface_real = ((medians['beta_a']+s29z1['redshift_BEAGLE']*medians['beta_b'])*s29z1['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+s29z1['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+s29z1['redshift_BEAGLE']*medians['beta_b']))

log_p_xi_eta_theta = norm.logpdf(s29z1['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(s29z1['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
p_bad = norm.pdf(s29z1['sfr_BEAGLE_instant'], scale=medians['outlier_sigma'], loc=medians['outlier_mean'])

z_bad = medians['pbad']*p_bad
z_good = (1-medians['pbad'])*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]


ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

scatter = ax1.scatter(s29z1['mass_BEAGLE_stellar'][idx_sort], s29z1['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5, s=100)

print(len(s29z1['mass_BEAGLE_stellar']))

xlow = 6.5
xhigh = 10.5
ylow = -2.5
yhigh = 2.5

ximin = 8.5
ximax = 10.0

lw = 3

ax1.plot((9.7,9.7), (ylow, yhigh), color='gray', linestyle='dashed', linewidth=2)

x_tmp = np.array([xlow, xhigh])
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, 1.25 $<$ z $<$ 2.0')
 
sig = rc_c1k['sig0_50_arr'][0] * ( ((1.0-rc_c1k['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)

# SANTINI, log(SFR) = alpha log(M / M_9p7) + beta
# 1.3 < z < 2
s1_alpha = 1.04
s1_beta = 1.01
s1_alpha_err = 0.03
s1_beta_err = 0.04

# santini z=1
s_x = np.array([8.4, 9.2, 10.0])
s_y = s1_alpha*(s_x - 9.7) + s1_beta
s_y_intrinsic = np.array([0.36, 0.35, 0.0])
s_y_observed = np.array([0.51, 0.46, 0.26])

ax1.plot(s_x, s_y, color='k', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
ax1.plot(s_x, s_y+s_y_intrinsic, color='k', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
ax1.plot(s_x, s_y-s_y_intrinsic, color='k', linestyle='--', linewidth=lw)
ax1.plot(s_x, s_y+s_y_observed, color='k', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
ax1.plot(s_x, s_y-s_y_observed, color='k', linestyle=':', linewidth=lw)

ax1.set_xlim(xlow, xhigh)
ax1.set_xlabel(r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(ylow, yhigh)
ax1.set_ylabel(r'$\log(\mathrm{SFR} \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

ax1.legend(loc='lower right', frameon=True, fontsize=1.44*fontsize_legend, framealpha=1.0)


cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\mathrm{MS}\longleftarrow \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \longrightarrow \mathrm{Outliers}$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=0, width=2, direction='in', labelsize=0*fontsize_axes)




#save = False
if save:
    plt.savefig('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/EAS_poster/9.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()






#%%
# =============================================================================
# Santini+17 'True' values - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

# converting normalisation
alpha_san = B_san - 9.7*A_san
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_san = A_san
beta_err_san = A_err_san
alpha_san_n = alpha_san + (normalisation*beta_san) # santini normalised
alpha_err_san_n = (alpha_err_san**2 - (normalisation*beta_err_san)**2) ** 0.5

# =============================================================================
# Santini+17 Original values - obtained by eye - delayed SFH, SFR from UV slope
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_san0 = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san0 = np.array((1.05, 1.1, 0.9, 0.75, 0.55))
A_err_san0 = np.array((0.03, 0.03, 0.04, 0.05, 0.18))
B_san0 = np.array((1.0, 1.15, 1.25, 1.2, 1.5))
B_err_san0 = np.array((0.05, 0.03, 0.03, 0.06, 0.12))

# converting normalisation
alpha_san0 = B_san0 - 9.7*A_san0
alpha_err_san0 = (B_err_san0**2 + (9.7*A_err_san0)**2) ** 0.5
beta_san0 = A_san0
beta_err_san0 = A_err_san0
alpha_san0_n = alpha_san0 + (normalisation*beta_san0) # santini normalised
alpha_err_san0_n = (alpha_err_san0**2 - (normalisation*beta_err_san0)**2) ** 0.5

# =============================================================================
# Speagle+14 - errors calculated dodgily
# =============================================================================
# log SFR(M∗, t) = (0.84 ± 0.02 − 0.026 ± 0.003 × t ) logM∗−(6.51 ± 0.24 − 0.11 ± 0.03 × t ), where t is the age of the universe in Gyr.

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z_speagle = np.linspace(0.5, 6.5, 1000)
t_speagle = cosmo.age(z_speagle).value # Gyr

alpha_speagle = -(6.51 - (0.11*t_speagle) )
alpha_err_speagle = 0.24 + (0.03*t_speagle)
beta_speagle = 0.84 - (0.026*t_speagle)
beta_err_speagle = 0.02 + (0.003*t_speagle)
alpha_speagle_n = alpha_speagle + (normalisation*beta_speagle) # santini normalised

# =============================================================================
# Schreiber+15 - ignoring high mass
# =============================================================================
#r ≡ log10(1 + z) and m ≡ log10(M∗/10**9 M):
#log10(SFRMS[M/yr]) = m − m0 + a0r − a1 max(0,m − m1 − a2 r)**2, 
#with m0 = 0.5 ± 0.07, a0 = 1.5 ± 0.15, a1 = 0.3 ± 0.08, m1 = 0.36 ± 0.3 and a2 = 2.5 ± 0.6.

z_schreiber = np.linspace(0.5, 4.0, 1000)
r_schreiber = np.log10(1+z_schreiber)

m0_schreiber = 0.5
a0_schreiber = 1.5
#a1 = 0.3
#m1 = 0.36
#a2 = 2.5

# m - m1 - a2r is usually < 0, except high mass, low redshift, IGNORED FOR NOW
#print( np.log10((10**9.9)/(1e9))- m1 - (a2*r_schreiber))

alpha_schreiber = - (9.0 + m0_schreiber - (a0_schreiber*r_schreiber))
beta_schreiber = np.linspace(1.0, 1.0, 1000)
alpha_schreiber_n = alpha_schreiber + (normalisation*beta_schreiber) # santini normalised

# =============================================================================
# Salmon+15
# =============================================================================
#log(SFR/M yr−1) = a log(M/M) + b
z_salmon = np.array((4.0, 5.0, 6.0))

alpha_salmon = np.array((-5.7, -4.4, -3.9))
alpha_err_salmon = np.array((2.1, 2.6, 1.6))
beta_salmon = np.array((0.7, 0.59, 0.54))
beta_err_salmon = np.array((0.21, 0.26, 0.16))
alpha_salmon_n = alpha_salmon + (normalisation*beta_salmon) # santini normalised

# =============================================================================
# Steinhardt+14
# =============================================================================
#log SFR(M yr−1) = α × (logM∗/M − 10.5) + β,
z_steinhardt = np.array(((4.0 + 4.8)/2, (4.8 + 6.0)/2))

beta_steinhardt = np.array((0.78, 0.78))
beta_err_steinhardt = np.array((0.02, 0.02))
alpha_steinhardt = np.array((1.976, 2.110)) - (10.5*beta_steinhardt)
alpha_steinhardt_n = alpha_steinhardt + (normalisation*beta_steinhardt) # santini normalised

# =============================================================================
# Tomczak+16 - Not obvious how to get these, also maybe not just star forming?!
# =============================================================================

# =============================================================================
# Schreiber+16 - single point, hard to find paper
# =============================================================================

# =============================================================================
# Kurczynski+16 - nice table, also mentions purpose of rescaling which I haven't included here
# =============================================================================
#log SFR = a ´ log M* + b + N (0, sIS).
z_kurc = np.array(((0.5 + 1.0)/2, (1.0 + 1.5)/2, (1.5 + 2.0)/2, (2.0 + 2.5)/2, (2.5 + 3.0)/2))

alpha_kurc = np.array((-8.394, -7.474, -7.484, -7.513, -7.729))
alpha_err_kurc = np.array((0.011, 0.010, 0.011, 0.018, 0.015))
beta_kurc = np.array((0.919, 0.825, 0.867, 0.849, 0.899))
beta_err_kurc = np.array((0.017, 0.012, 0.013, 0.021, 0.017))
alpha_kurc_n = alpha_kurc + (normalisation*beta_kurc) # santini normalised

# =============================================================================
# THE PLOT
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = 0.3
xhigh = 6.7

param = ['beta', 'alpha', 'sig0']

ax1 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax3 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1,ax2,ax3]

ax1.scatter(z_salmon, beta_salmon, label='Salmon+15')
ax1.errorbar(z_salmon, beta_salmon, yerr=beta_err_salmon, ls='none')
ax1.scatter(z_steinhardt, beta_steinhardt, label='Steinhardt+14')
ax1.errorbar(z_steinhardt, beta_steinhardt, yerr=beta_err_steinhardt, ls='none')
ax1.scatter(z_san0, beta_san0, label='Santini+17 Observed')
ax1.errorbar(z_san0, beta_san0, yerr=beta_err_san0, ls='none')
ax1.scatter(z_san, beta_san, label='Santini+17 Simulation')
ax1.errorbar(z_san, beta_san, yerr=beta_err_san, ls='none')
ax1.scatter(z_kurc, beta_kurc, label='Kurczynski+16')
ax1.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none')
ax1.plot(z_speagle, beta_speagle, label='Speagle+14', linestyle=':', color='gray')
ax1.plot(z_schreiber, beta_schreiber, label='Schreiber+15', linestyle='dashed', color='gray')

alpha_err_salmon_n = np.zeros(len(alpha_salmon_n))
alpha_err_steinhardt_n = np.zeros(len(alpha_steinhardt_n))
alpha_err_kurc_n = np.zeros(len(alpha_kurc_n))


ax2.scatter(z_salmon, alpha_salmon_n, label='Salmon+15')
ax2.errorbar(z_salmon, alpha_salmon_n, yerr=alpha_err_salmon_n, ls='none')
ax2.scatter(z_steinhardt, alpha_steinhardt_n, label='Steinhardt+14')
ax2.errorbar(z_steinhardt, alpha_steinhardt_n, yerr=alpha_err_steinhardt_n, ls='none')
ax2.scatter(z_san0, alpha_san0_n, label='Santini+17 Observed')
ax2.errorbar(z_san0, alpha_san0_n, yerr=alpha_err_san0_n, ls='none')
ax2.scatter(z_san, alpha_san_n, label='Santini+17 Simulation')
ax2.errorbar(z_san, alpha_san_n, yerr=alpha_err_san_n, ls='none')
ax2.scatter(z_kurc, alpha_kurc_n, label='Kurczynski+16')
ax2.errorbar(z_kurc, alpha_kurc_n, yerr=alpha_err_kurc_n, ls='none')
ax2.plot(z_speagle, alpha_speagle_n, label='Speagle+14', linestyle=':', color='gray')
ax2.plot(z_schreiber, alpha_schreiber_n, label='Schreiber+15', linestyle='dashed', color='gray')

ax3.scatter(0, 0, label='Salmon+15')
ax3.scatter(0, 0, label='Steinhardt+14')
ax3.scatter(0, 0, label='Santini+17 Observed')
ax3.scatter(0, 0, label='Santini+17 Simulation')
ax3.scatter(0, 0, label='Kurczynski+16')
ax3.plot(0, 0, label='Our Work, 1.25 $<$ z $<$ 6.0', color='red')
ax3.plot(0, 0, label='Our Work, redshift bins', color='blue')
ax3.plot(0, 0, label='Speagle+14', linestyle=':', color='gray')
ax3.plot(0, 0, label='Schreiber+15', linestyle='dashed', color='gray')



for ax in [ax1,ax2,ax3]:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax.xaxis.set_tick_params(which='minor', size=0, width=2, direction='in', bottom='on', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(which='minor', size=0, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

# redshift bins
#for rc in [rc_268, rc_269, rc_270, rc_273]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k')
##        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='blue', zorder=0)

for rc in [rc_c1, rc_c2, rc_c3, rc_c4, rc_c5]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='blue')
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.1, color='blue', zorder=0)
        
for rc in [rc_c]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='red')
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='red', zorder=0)
        
#        
#for rc in [rc_cla]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', linestyle='dashed')
#        
#for rc in [rc_cp1, rc_cp2, rc_cp3, rc_cp4, rc_cp]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta')
#
#for rc in [rc_cpla]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', linestyle='dashed')        
#
#for rc in [rc_265, rc_274]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])])        

# layout
ylim_low = [0.2, 0.1, 0.1]
ylim_high = [1.4, 2.4, 0.9]
string = [string_slope, string_normalisation, string_scatter]
ticker_maj = [0.4, 0.5, 0.2]
ticker_min = [0.1, 0.2, 0.1]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k', label='original z bins')
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', label='z bins and full run (c)')
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', linestyle='dashed', label='full run linear alpha (c)')
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', label='z bins and full run (cp)')
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', linestyle='dashed', label='full run linear alpha (cp)')   
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='#1f77b4', label='original beta == 1')
#ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='#ff7f0e', label='original beta == const, z$<$4')
ax3.legend(loc='upper left', frameon=False, fontsize=1.22*fontsize_legend, ncol=2)


if save:
    plt.savefig('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/EAS_poster/4.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#%%








