#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 08:41:18 2021

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
fontsize_legend = 12
fontsize_axes = 20
figuresize = 7

normalisation = 9.7

load = True
#load = False
save = False
 
# =============================================================================
# useful axis strings
# =============================================================================

string_slope = r'$\mathrm{MS\,Slope,}\,\beta$'
string_normalisation = r'MS Normalisation, $\alpha$'
string_scatter = r'$\mathrm{MS\,Scatter,}\,\sigma$'
string_ssfr = r'$\log(\mathrm{sSFR}/\mathrm{yr}^{-1})$'
string_mass = r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$'
string_sfr = r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$'
string_deltaMS = r'$\Delta_{MS}$'
string_prob_ratio = r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$'
string_bias_test = r'$\Delta \mathrm{Parameter}$'
string_pbad = r'pbad'
string_outlier_mean = r'outlier mean'
string_outlier_sigma = r'outlier sigma'

# =============================================================================
# opening data
# =============================================================================
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3

if load:


    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
        chain_MS_29_c1i = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
        chain_MS_29_c2i = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
        chain_MS_29_c3i = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z4p0-5p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
        chain_MS_29_c4i = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z5p0-6p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all.p', 'rb') as f:
        chain_MS_29_c5i = pickle.load(f, encoding='latin1') 


    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x5000.p', 'rb') as f:
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
        
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z1p25-2p0_4x5000.p', 'rb') as f:
        chain_MS_29_cp1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z2p0-3p0_4x5000.p', 'rb') as f:
        chain_MS_29_cp2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z3p0-4p0_4x5000.p', 'rb') as f:
        chain_MS_29_cp3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z4p0-5p0_4x5000.p', 'rb') as f:
        chain_MS_29_cp4 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z5p0-6p0_4x5000.p', 'rb') as f:
        chain_MS_29_cp5 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters+parallels_z1p25-6p0_4x5000_ssfr_alpha.p', 'rb') as f:
        chain_MS_29_cp = pickle.load(f, encoding='latin1')

    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z1-2.p', 'rb') as f:
        chain_MS_29_c1f = pickle.load(f, encoding='latin1')        
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z2-3.p', 'rb') as f:
        chain_MS_29_c2f = pickle.load(f, encoding='latin1')             
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z3-4.p', 'rb') as f:
        chain_MS_29_c3f = pickle.load(f, encoding='latin1')             
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_z4-5.p', 'rb') as f:
        chain_MS_29_c4f = pickle.load(f, encoding='latin1')             

    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-6p0_4x5000_ssfr_alpha_fixed_pbad.p', 'rb') as f:
        chain_MS_29_c_pbad_bins = pickle.load(f, encoding='latin1')        
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z3p0-6p0_4x5000_z4-5.p', 'rb') as f:
        chain_MS_29_c4h = pickle.load(f, encoding='latin1')
        
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_k_fitted.p', 'rb') as f:
        chain_MS_29_c1_5000k = pickle.load(f, encoding='latin1')          
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x5000.p', 'rb') as f:
        chain_MS_29_c1_5000 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_k_fitted_1_2.p', 'rb') as f:
        chain_MS_29_c1_10000k12 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_k_fitted.p', 'rb') as f:
        chain_MS_29_c1_10000k = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x10000.p', 'rb') as f:
        chain_MS_29_c1_10000 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_true_constant_alpha.p', 'rb') as f:
        chain_MS_29_c1_50000t = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_29_clusters_z1p25-2p0_4x50000.p', 'rb') as f:
        chain_MS_29_c1_50000 = pickle.load(f, encoding='latin1')   
        




        
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z1p25-2p0.p', 'rb') as f:
        chain_MS_27_c1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z2p0-3p0.p', 'rb') as f:
        chain_MS_27_c2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z3p0-4p0.p', 'rb') as f:
        chain_MS_27_c3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z4p0-5p0.p', 'rb') as f:
        chain_MS_27_c4 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z1p25-6p0.p', 'rb') as f:
        chain_MS_27_c = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters_z1p25-6p0_linear_alphaN.p', 'rb') as f:
        chain_MS_27_cla = pickle.load(f, encoding='latin1')
        
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z1p25-2p0.p', 'rb') as f:
        chain_MS_27_cp1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z2p0-3p0.p', 'rb') as f:
        chain_MS_27_cp2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z3p0-4p0.p', 'rb') as f:
        chain_MS_27_cp3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z4p0-5p0.p', 'rb') as f:
        chain_MS_27_cp4 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z1p25-6p0.p', 'rb') as f:
        chain_MS_27_cp = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s27_clusters+parallels_z1p25-6p0_linear_alphaN.p', 'rb') as f:
        chain_MS_27_cpla = pickle.load(f, encoding='latin1')    
    
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_268_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_268 = pickle.load(f, encoding='latin1') #2.0
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_269_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_269 = pickle.load(f, encoding='latin1') #3.0
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_270_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_270 = pickle.load(f, encoding='latin1') #4.0
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_265_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_265 = pickle.load(f, encoding='latin1') #beta=1
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_274_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_274 = pickle.load(f, encoding='latin1') #beta=constant from z<4
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_273_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_273 = pickle.load(f, encoding='latin1') #1.25    
    
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_27_clusters_z1p25-6p0_4x5000_test4_z4_run1.p', 'rb') as f:
        chain_MS_27_c4_test4_run1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_27_clusters_z1p25-6p0_4x2000_test4_z4_run2.p', 'rb') as f:
        chain_MS_27_c4_test4_run2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_scenario_27_clusters_z1p25-6p0_4x2000_test4_z4_run4.p', 'rb') as f:
        chain_MS_27_c4_test4_run4 = pickle.load(f, encoding='latin1')  
        
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


        
if False:

    # converged chains
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_154_subset_002.p', 'rb') as f:
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_265_subset_001.p', 'rb') as f: # 9000-10000 
#        chain_MS = pickle.load(f, encoding='latin1') 
#    print(chain_burn.dtype.names)

    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_267_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_267 = pickle.load(f, encoding='latin1') #1.0

    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_278_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_278 = pickle.load(f, encoding='latin1') #k!=1

    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z1p25-2p0.p', 'rb') as f:
        chain_MS_26_c1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z2p0-3p0.p', 'rb') as f:
        chain_MS_26_c2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z3p0-4p0.p', 'rb') as f:
        chain_MS_26_c3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z4p0-5p0.p', 'rb') as f:
        chain_MS_26_c4 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z5p0-6p0.p', 'rb') as f:
        chain_MS_26_c5 = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z1p25-6p0.p', 'rb') as f:
        chain_MS_26_c = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z1p25-6p0_z_upper_5.p', 'rb') as f:
        chain_MS_26_cu5 = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters_z1p25-6p0_z_upper_4.p', 'rb') as f:
        chain_MS_26_cu4 = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z1p25-2p0.p', 'rb') as f:
        chain_MS_26_cp1 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z2p0-3p0.p', 'rb') as f:
        chain_MS_26_cp2 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z3p0-4p0.p', 'rb') as f:
        chain_MS_26_cp3 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z4p0-5p0.p', 'rb') as f:
        chain_MS_26_cp4 = pickle.load(f, encoding='latin1')  
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z5p0-6p0.p', 'rb') as f:
        chain_MS_26_cp5 = pickle.load(f, encoding='latin1')   
    with open('/Users/lester/Documents/linmix_files/lm_chain_s26_clusters+parallels_z1p25-6p0.p', 'rb') as f:
        chain_MS_26_cp = pickle.load(f, encoding='latin1')

#    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p', 'rb') as f:
#        astrodeep_pickle = pickle.load(f, encoding='latin1') 
#    print(astrodeep_pickle.keys())



    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
        data = pickle.load(f, encoding='latin1') 
#    print(data.keys())

    bias_tests = ['110', '111', '112', '113', '114', '115', '116', '117', '118', '119']
    bias_tests = ['115', '117', '110', '119', '112', '111', '116', '113', '118', '114']
    chain_bias = []
    for test in bias_tests:
        with open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_truncated_{}.p'.format(test), 'rb') as f:
            chain_bias.append(pickle.load(f, encoding='latin1'))



    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.fits'
    s23 = fits.open(fileName)
    print(s23.info())
    print(s23[1].header)
    s23 = fits.open(fileName)[1].data # scenario 23


#%%
#def read_chain(zlow, zhigh, chain_MS, fit):
#    
#    z = np.linspace(zlow, zhigh, 1000)
#
#    redshift_arr = np.repeat(np.array([z]).T, len(chain_MS), axis=1).T
#
#    beta_a_arr = np.array([chain_MS['beta_a']]).T
#    beta_b_arr = np.array([chain_MS['beta_b']]).T
#    beta_arr = beta_a_arr + redshift_arr*beta_b_arr
#    beta_16_arr = np.percentile(beta_arr, 16, axis=0)
#    beta_50_arr = np.median(beta_arr, axis=0)
#    beta_84_arr = np.percentile(beta_arr, 84, axis=0)
#
#    if fit == 'const_beta_ssfr_alpha' or fit == 'const_beta_ssfr_alpha_fixed_pbad':
#
#        ssfr_a_arr = np.array([chain_MS['alphaN_a']]).T
#        ssfr_b_arr = np.array([chain_MS['alphaN_b']]).T
#        ssfr_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr) - 9.0
#        ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
#        ssfr_50_arr = np.median(ssfr_arr, axis=0)
#        ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)
#        
#        alpha_arr = ssfr_arr + normalisation
#        alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
#        alpha_50_arr = np.median(alpha_arr, axis=0)
#        alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)
#        
#    elif fit == 'const_beta_linear_alpha':
#
#        alpha_a_arr = np.array([chain_MS['alphaN_a']]).T
#        alpha_b_arr = np.array([chain_MS['alphaN_b']]).T        
#        alpha_arr = alpha_a_arr + redshift_arr*alpha_b_arr
#        alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
#        alpha_50_arr = np.median(alpha_arr, axis=0)
#        alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)
#
#        ssfr_arr = alpha_arr - normalisation
#        ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
#        ssfr_50_arr = np.median(ssfr_arr, axis=0)
#        ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)
#        
#    if fit == 'const_beta_ssfr_alpha_fixed_pbad':
#        pbad_arr = np.array([chain_MS['pbad'][:,0]]).T + redshift_arr*0
#        pbad_16_arr = np.percentile(pbad_arr, 16, axis=0)
#        pbad_50_arr = np.median(pbad_arr, axis=0)
#        pbad_84_arr = np.percentile(pbad_arr, 84, axis=0)
#    else:
#        pbad_arr = np.array([chain_MS['pbad']]).T + redshift_arr*0
#        pbad_16_arr = np.percentile(pbad_arr, 16, axis=0)
#        pbad_50_arr = np.median(pbad_arr, axis=0)
#        pbad_84_arr = np.percentile(pbad_arr, 84, axis=0)
#        
#        
##        print(chain_MS['pbad'][:,0])
##        chain_MS['pbad'] = chain_MS['pbad'][:,0]
##        exit
#        
#    sig0_arr = np.array([chain_MS['sig0']]).T + redshift_arr*0
#    sig0_16_arr = np.percentile(sig0_arr, 16, axis=0)
#    sig0_50_arr = np.median(sig0_arr, axis=0)
#    sig0_84_arr = np.percentile(sig0_arr, 84, axis=0)
#
#    outlier_mean_arr = np.array([chain_MS['outlier_mean']]).T + redshift_arr*0
#    outlier_mean_16_arr = np.percentile(outlier_mean_arr, 16, axis=0)
#    outlier_mean_50_arr = np.median(outlier_mean_arr, axis=0)
#    outlier_mean_84_arr = np.percentile(outlier_mean_arr, 84, axis=0)
#    
#    outlier_sigma_arr = np.array([chain_MS['outlier_sigma']]).T + redshift_arr*0
#    outlier_sigma_16_arr = np.percentile(outlier_sigma_arr, 16, axis=0)
#    outlier_sigma_50_arr = np.median(outlier_sigma_arr, axis=0)
#    outlier_sigma_84_arr = np.percentile(outlier_sigma_arr, 84, axis=0)
#    
#    # DICTIONARY
#    
#    dic = {}
#    
#    dic['z'] = z # array of 1000 numbers between zlow and zhigh
#    
#    dic['beta_16_arr'] = beta_16_arr # 1000 values (1 per redshift interval)
#    dic['beta_50_arr'] = beta_50_arr # 1000 values (1 per redshift interval)
#    dic['beta_84_arr'] = beta_84_arr # 1000 values (1 per redshift interval)
#    
#    dic['ssfr_16_arr'] = ssfr_16_arr # 1000 values (1 per redshift interval)
#    dic['ssfr_50_arr'] = ssfr_50_arr # 1000 values (1 per redshift interval)
#    dic['ssfr_84_arr'] = ssfr_84_arr # 1000 values (1 per redshift interval)
#
#    dic['alpha_16_arr'] = alpha_16_arr # 1000 values (1 per redshift interval)
#    dic['alpha_50_arr'] = alpha_50_arr # 1000 values (1 per redshift interval)
#    dic['alpha_84_arr'] = alpha_84_arr # 1000 values (1 per redshift interval)
#    
#    dic['sig0_16_arr'] = sig0_16_arr # 1000 values (1 per redshift interval)
#    dic['sig0_50_arr'] = sig0_50_arr # 1000 values (1 per redshift interval)
#    dic['sig0_84_arr'] = sig0_84_arr # 1000 values (1 per redshift interval)
#
#    dic['pbad_16_arr'] = pbad_16_arr # 1000 values (1 per redshift interval)
#    dic['pbad_50_arr'] = pbad_50_arr # 1000 values (1 per redshift interval)
#    dic['pbad_84_arr'] = pbad_84_arr # 1000 values (1 per redshift interval)
#
#    dic['outlier_mean_16_arr'] = outlier_mean_16_arr # 1000 values (1 per redshift interval)
#    dic['outlier_mean_50_arr'] = outlier_mean_50_arr # 1000 values (1 per redshift interval)
#    dic['outlier_mean_84_arr'] = outlier_mean_84_arr # 1000 values (1 per redshift interval)    
#    
#    dic['outlier_sigma_16_arr'] = outlier_sigma_16_arr # 1000 values (1 per redshift interval)
#    dic['outlier_sigma_50_arr'] = outlier_sigma_50_arr # 1000 values (1 per redshift interval)
#    dic['outlier_sigma_84_arr'] = outlier_sigma_84_arr # 1000 values (1 per redshift interval)    
#    
#    return dic

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

rc_268 = read_chain(2.0, 3.0, chain_MS_268, 'const_beta_ssfr_alpha') #2.0
rc_269 = read_chain(3.0, 4.0, chain_MS_269, 'const_beta_ssfr_alpha')
rc_270 = read_chain(4.0, 5.0, chain_MS_270, 'const_beta_ssfr_alpha')
rc_273 = read_chain(1.25, 2.0, chain_MS_273, 'const_beta_ssfr_alpha')

rc_265 = read_chain(1.25, 6.0, chain_MS_265, 'const_beta_ssfr_alpha') # beta=1
rc_274 = read_chain(1.25, 6.0, chain_MS_274, 'const_beta_ssfr_alpha') # beta=constant from z<4

rc_c1 = read_chain(1.25, 2.0, chain_MS_29_c1, 'const_beta_ssfr_alpha')
rc_c2 = read_chain(2.0, 3.0, chain_MS_29_c2, 'const_beta_ssfr_alpha')
rc_c3 = read_chain(3.0, 4.0, chain_MS_29_c3, 'const_beta_ssfr_alpha')
rc_c4 = read_chain(4.0, 5.0, chain_MS_29_c4, 'const_beta_ssfr_alpha')
rc_c5 = read_chain(5.0, 6.0, chain_MS_29_c5, 'const_beta_ssfr_alpha')
rc_c = read_chain(1.25, 6.0, chain_MS_29_c, 'const_beta_ssfr_alpha')


rc_c1_5000k = read_chain(1.25, 2.0, chain_MS_29_c1_5000k, 'const_beta_ssfr_alpha')
rc_c1_5000 = read_chain(1.25, 2.0, chain_MS_29_c1_5000, 'const_beta_ssfr_alpha')
rc_c1_10000k12 = read_chain(1.25, 2.0, chain_MS_29_c1_10000k12, 'const_beta_ssfr_alpha')
rc_c1_10000k = read_chain(1.25, 2.0, chain_MS_29_c1_10000k, 'const_beta_ssfr_alpha')
rc_c1_10000 = read_chain(1.25, 2.0, chain_MS_29_c1_10000, 'const_beta_ssfr_alpha')
rc_c1_50000t = read_chain(1.25, 2.0, chain_MS_29_c1_50000t, 'const_beta_linear_alpha') # const_beta_linear_alpha
rc_c1_50000 = read_chain(1.25, 2.0, chain_MS_29_c1_50000, 'const_beta_ssfr_alpha')


rc_cp1 = read_chain(1.25, 2.0, chain_MS_27_cp1, 'const_beta_ssfr_alpha')
rc_cp2 = read_chain(2.0, 3.0, chain_MS_27_cp2, 'const_beta_ssfr_alpha')
rc_cp3 = read_chain(3.0, 4.0, chain_MS_27_cp3, 'const_beta_ssfr_alpha')
rc_cp4 = read_chain(4.0, 5.0, chain_MS_27_cp4, 'const_beta_ssfr_alpha')
rc_cp5 = read_chain(5.0, 6.0, chain_MS_29_cp5, 'const_beta_ssfr_alpha')
rc_cp = read_chain(1.25, 6.0, chain_MS_27_cp, 'const_beta_ssfr_alpha')


rc_c1f = read_chain(1.25, 2.0, chain_MS_29_c1f, 'const_beta_ssfr_alpha')
rc_c2f = read_chain(2.0, 3.0, chain_MS_29_c2f, 'const_beta_ssfr_alpha')
rc_c3f = read_chain(3.0, 4.0, chain_MS_29_c3f, 'const_beta_ssfr_alpha')
rc_c4f = read_chain(4.0, 5.0, chain_MS_29_c4f, 'const_beta_ssfr_alpha')
rc_c4h = read_chain(4.0, 5.0, chain_MS_29_c4h, 'const_beta_ssfr_alpha')

rc_c_pbad_bins = read_chain(1.25, 6.0, chain_MS_29_c_pbad_bins, 'const_beta_ssfr_alpha_fixed_pbad')


rc_c4_test4_run1 = read_chain(4.0, 5.0, chain_MS_27_c4_test4_run1, 'const_beta_ssfr_alpha')
rc_c4_test4_run2 = read_chain(4.0, 5.0, chain_MS_27_c4_test4_run2, 'const_beta_ssfr_alpha')
rc_c4_test4_run4 = read_chain(4.0, 5.0, chain_MS_27_c4_test4_run4, 'const_beta_ssfr_alpha')

rc_c1i = read_chain(1.25, 2.0, chain_MS_29_c1i, 'const_beta_linear_alpha')
rc_c2i = read_chain(2.0, 3.0, chain_MS_29_c2i, 'const_beta_linear_alpha')
rc_c3i = read_chain(3.0, 4.0, chain_MS_29_c3i, 'const_beta_linear_alpha')
rc_c4i = read_chain(4.0, 5.0, chain_MS_29_c4i, 'const_beta_linear_alpha')
rc_c5i = read_chain(5.0, 6.0, chain_MS_29_c5i, 'const_beta_linear_alpha')


#print(rc_c3['alpha_50_arr'][0])
#print(rc_c3f['alpha_50_arr'][0])
#print(rc_c4['alpha_50_arr'][0])
#print(rc_c4f['alpha_50_arr'][0])
#
#print(rc_c3['beta_50_arr'][0])
#print(rc_c3f['beta_50_arr'][0])
#print(rc_c4['beta_50_arr'][0])
#print(rc_c4f['beta_50_arr'][0])
#
#print(rc_c4h['alpha_50_arr'][0])
#print(rc_c4h['beta_50_arr'][0])

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

medians = get_medians(chain_MS_27_c)



#%%

# =============================================================================
# Parameter vs redshift, z1-2 tests
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = 0.7
xhigh = 9.3

param = ['beta', 'alpha', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']

ax1 = fig.add_axes([0, 2.5, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 2.0, 0.5, 0.5]) #[left, bottom, width, height]
ax3 = fig.add_axes([0, 1.5, 0.5, 0.5]) #[left, bottom, width, height]
ax4 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]
ax5 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax6 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1,ax2,ax3,ax4,ax5,ax6]

for ax in axes:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

# redshift bins
count = 0
#for rc in [rc_c1_5000k, rc_c1_5000, rc_c1_10000k12, rc_c1_10000k, rc_c1_10000, rc_c1_50000t, rc_c1_50000]:
for rc in [rc_c1_5000,rc_c1_10000, rc_c1_50000, rc_c1_50000t, rc_c1_5000k, rc_c1_10000k, rc_c1_10000k12, rc_c1i]:
#for rc in [rc_c1_5000k, rc_c1_5000, rc_c1_50000]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z']+count, rc['{}_50_arr'.format(param[a])])
        ax.fill_between(rc['z']+count, rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, zorder=0)
    count+=1

for a, ax in enumerate(axes):
    ax.plot(0, 0, color='#1f77b4', label='5000')
    ax.plot(0, 0, color='#ff7f0e', label='10000')
    ax.plot(0, 0, color='#2ca02c', label='50000')
    ax.plot(0, 0, color='#d62728', label='50000t')
    ax.plot(0, 0, color='#9467bd', label='5000k')
    ax.plot(0, 0, color='#8c564b', label='10000k')
    ax.plot(0, 0, color='#e377c2', label='10000k12')
    ax.plot(0, 0, color='gray', label='new')
    ax.legend(loc='upper left', frameon=False, fontsize=fontsize_legend)

# layout
ylim_low = [0.6, 0.7, 0.2, 0.0, 0.6, 0.6]
ylim_high = [1.0, 1.1, 0.4, 0.3, 1.1, 1.1]
string = [string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma]
ticker_maj = [0.04, 0.1, 0.04, 0.04, 0.1, 0.04]
ticker_min = [0.04, 0.1, 0.04, 0.04, 0.1, 0.04]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

#ax6.set_xlabel('Redshift', labelpad=10)
ax6.xaxis.set_tick_params(labelsize=0*fontsize_axes)

#ax6.plot(0, 0, color='blue', label='redshift bins')
#ax6.plot(0, 0, color='red', label='full run')
#ax6.plot(0, 0, color='k', linestyle='dashed', label='original z bins')
#ax6.plot(0, 0, color='orange', label='full sample')
#ax6.plot(0, 0, color='orange', linestyle='dashed', label='z3-6 sample')
#ax6.plot(0, 0, color='k', label='pbad fixed from z bins')
#ax6.legend(bbox_to_anchor=(0.0, 1.0), loc=2, frameon=False, fontsize=fontsize_legend)

if save:
    plt.savefig('014_redshift_bins.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%



#%%

# =============================================================================
# Parameter vs redshift, in redshift bins with new method 
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = 0.3
xhigh = 6.7

param = ['beta', 'alpha', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']

ax1 = fig.add_axes([0, 2.5, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 2.0, 0.5, 0.5]) #[left, bottom, width, height]
ax3 = fig.add_axes([0, 1.5, 0.5, 0.5]) #[left, bottom, width, height]
ax4 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]
ax5 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax6 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1,ax2,ax3,ax4,ax5,ax6]

for ax in axes:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

# redshift bins
for rc in [rc_268, rc_269, rc_270, rc_273]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k', linestyle='dashed')
#        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='blue', zorder=0)

for rc in [rc_c1i, rc_c2i, rc_c3i, rc_c4i, rc_c5i]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='blue')
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='blue', zorder=0)

for rc in [rc_c]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='red')
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='red', zorder=0)

#TEST 4
#for rc in [rc_c4_test4_run1, rc_c4_test4_run2, rc_c4_test4_run4]:
#    for a, ax in enumerate(axes):
#        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='orange')
#        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='orange', zorder=0)

#TEST 4
for rc in [rc_c1f, rc_c2f, rc_c3f, rc_c4f]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='orange')
#        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='orange', zorder=0)
for rc in [rc_c4h]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='orange', linestyle='dashed')
#        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='orange', zorder=0)
        
#Full Run but fixing pbad to the value of the redshift bin results
for rc in [rc_c_pbad_bins]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k')
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='k', zorder=0)


# layout
ylim_low = [0.5, 0.2, 0.05, -0.05, -1.5, 0.0]
ylim_high = [1.3, 2.2, 0.45, 0.45, 3.5, 5.5]
string = [string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma]
ticker_maj = [0.2, 0.4, 0.1, 0.1, 1.0, 1.0]
ticker_min = [0.1, 0.2, 0.05, 0.05, 0.5, 0.5]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

ax6.set_xlabel('Redshift', labelpad=10)
ax6.xaxis.set_tick_params(labelsize=fontsize_axes)

ax6.plot(0, 0, color='blue', label='redshift bins')
ax6.plot(0, 0, color='red', label='full run')
ax6.plot(0, 0, color='k', linestyle='dashed', label='original z bins')
ax6.plot(0, 0, color='orange', label='full sample')
ax6.plot(0, 0, color='orange', linestyle='dashed', label='z3-6 sample')
ax6.plot(0, 0, color='k', label='pbad fixed from z bins')
ax6.legend(bbox_to_anchor=(0.0, 1.0), loc=2, frameon=False, fontsize=fontsize_legend)

if save:
    plt.savefig('014_redshift_bins.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%


plt.plot(chain_MS_29_c_pbad_bins['pbad'][:,2])
plt.show()
plt.plot(chain_MS_29_c_pbad_bins['zeta'][:,2])
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

param = ['beta', 'alpha', 'ssfr']

ax1 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax3 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1,ax2,ax3]

'''
ax1.scatter(z_san, beta_san, label='Santini+17')
ax1.errorbar(z_san, beta_san, yerr=beta_err_san, ls='none')
ax1.scatter(z_san0, beta_san0, label='Santini+17 Raw')
ax1.errorbar(z_san0, beta_san0, yerr=beta_err_san0, ls='none')
ax1.scatter(z_salmon, beta_salmon, label='Salmon+15')
ax1.errorbar(z_salmon, beta_salmon, yerr=beta_err_salmon, ls='none')
ax1.scatter(z_steinhardt, beta_steinhardt, label='Steinhardt+14')
ax1.errorbar(z_steinhardt, beta_steinhardt, yerr=beta_err_steinhardt, ls='none')
ax1.scatter(z_kurc, beta_kurc, label='Kurczynski+16')
ax1.errorbar(z_kurc, beta_kurc, yerr=beta_err_kurc, ls='none')
ax1.plot(z_speagle, beta_speagle, label='Speagle+14', linestyle=':')
ax1.plot(z_schreiber, beta_schreiber, label='Schreiber+15', linestyle='dashed')
'''
'''
alpha_err_salmon_n = np.zeros(len(alpha_salmon_n))
alpha_err_steinhardt_n = np.zeros(len(alpha_steinhardt_n))
alpha_err_kurc_n = np.zeros(len(alpha_kurc_n))

ax2.scatter(z_san, alpha_san_n, label='Santini+17')
ax2.errorbar(z_san, alpha_san_n, yerr=alpha_err_san_n, ls='none')
ax2.scatter(z_san0, alpha_san0_n, label='Santini+17 Raw')
ax2.errorbar(z_san0, alpha_san0_n, yerr=alpha_err_san0_n, ls='none')
ax2.scatter(z_salmon, alpha_salmon_n, label='Salmon+15')
ax2.errorbar(z_salmon, alpha_salmon_n, yerr=alpha_err_salmon_n, ls='none')
ax2.scatter(z_steinhardt, alpha_steinhardt_n, label='Steinhardt+14')
ax2.errorbar(z_steinhardt, alpha_steinhardt_n, yerr=alpha_err_steinhardt_n, ls='none')
ax2.scatter(z_kurc, alpha_kurc_n, label='Kurczynski+16')
ax2.errorbar(z_kurc, alpha_kurc_n, yerr=alpha_err_kurc_n, ls='none')
ax2.plot(z_speagle, alpha_speagle_n, label='Speagle+14', linestyle=':')
ax2.plot(z_schreiber, alpha_schreiber_n, label='Schreiber+15', linestyle='dashed')
'''

for ax in [ax1,ax2,ax3]:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

# redshift bins
for rc in [rc_268, rc_269, rc_270, rc_273]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k')
#        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='blue', zorder=0)

for rc in [rc_c1, rc_c2, rc_c3, rc_c4, rc_c]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan')
        
for rc in [rc_cla]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', linestyle='dashed')
        
for rc in [rc_cp1, rc_cp2, rc_cp3, rc_cp4, rc_cp]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta')

for rc in [rc_cpla]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', linestyle='dashed')        

for rc in [rc_265, rc_274]:
    for a, ax in enumerate(axes):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])])        

# layout
ylim_low = [0.5, 0.2, -9.4]
ylim_high = [1.3, 2.2, -7.4]
string = [string_slope, string_normalisation, string_ssfr]
ticker_maj = [0.2, 0.4, 0.4]
ticker_min = [0.1, 0.2, 0.2]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='k', label='original z bins')
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', label='z bins and full run (c)')
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='cyan', linestyle='dashed', label='full run linear alpha (c)')
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', label='z bins and full run (cp)')
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='magenta', linestyle='dashed', label='full run linear alpha (cp)')   
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='#1f77b4', label='original beta == 1')
ax3.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='#ff7f0e', label='original beta == const, z$<$4')
ax3.legend(bbox_to_anchor=(0.0, 1.0), loc=2, frameon=False, fontsize=fontsize_legend)


if save:
    plt.savefig('001_parameters_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#%%

# =============================================================================
# new IRAC test (do 68 beagle fitted credible intervals overlap with +-1sigma measured input data)
# =============================================================================
ch1_beagle_flux_median = (10**6) * (10**((ADx['ch1_beagle_mag_median'] - 8.9)/(-2.5)))
ch1_beagle_flux_lower = (10**6) * (10**((ADx['ch1_beagle_mag_lower'] - 8.9)/(-2.5)))
ch1_beagle_flux_upper = (10**6) * (10**((ADx['ch1_beagle_mag_upper'] - 8.9)/(-2.5)))

ch2_beagle_flux_median = (10**6) * (10**((ADx['ch2_beagle_mag_median'] - 8.9)/(-2.5)))
ch2_beagle_flux_lower = (10**6) * (10**((ADx['ch2_beagle_mag_lower'] - 8.9)/(-2.5)))
ch2_beagle_flux_upper = (10**6) * (10**((ADx['ch2_beagle_mag_upper'] - 8.9)/(-2.5)))

b_CH1_AD = ADx['b_CH1_AD']
b_errCH1_AD = ADx['b_errCH1_AD']
b_CH2_AD = ADx['b_CH2_AD']
b_errCH2_AD = ADx['b_errCH2_AD']

# TRUE means REJECT because unreliable IRAC
# the checks are:
# 1sigma data DOES NOT overlap 68 credible interval
# can be the case in either CH1 OR CH2
# both CH1 & CH2 must also be -67 (unreliable from COVMAX > 1)
idx_new_IRAC = (((b_CH1_AD-b_errCH1_AD > ch1_beagle_flux_upper) | (b_CH1_AD+b_errCH1_AD < ch1_beagle_flux_lower)) | \
               ((b_CH2_AD-b_errCH2_AD > ch2_beagle_flux_upper) | (b_CH2_AD+b_errCH2_AD < ch2_beagle_flux_lower))) & \
               ((ADx['CH1_BEAGLE_input'] < -60) & (ADx['CH2_BEAGLE_input'] < -60))

# FALSE MEANS REJECT
# TRUE MEANS KEEP
idx_new_IRAC = ~idx_new_IRAC               


# =============================================================================
# redshift vs redshift inc IRAC flags
# =============================================================================

idx_clusters_parallels = 2.0 # clusters (for now)
z_lower = 1.25
z_upper = 6.0

idx1 = (ADx['field_AD']%idx_clusters_parallels==0.0) # 2 == 0 clusters, 2 == 1 parallels, 1 == 0 both
idx2 = (ADx['relflag_AD']==1.0) # relflag
idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

idx3_KIRAC = np.logical_and(ADx['CH1_BEAGLE_input']<-60.0, ADx['CH2_BEAGLE_input']<-60.0)
idx3_KIRAC = np.logical_and(idx3_KIRAC, ADx['Ks_BEAGLE_input']<-60.0)
idx3_KIRAC_removed = idx3_KIRAC
idx3_KIRAC = ~idx3_KIRAC

MCLmassLow = np.empty(ADx['redshift_BEAGLE'].shape)
for i in range(len(MCLmassLow)):
    if ADx['redshift_BEAGLE'][i] <2.1789654:
        MCLmassLow[i] = 8.0
    elif ADx['redshift_BEAGLE'][i] > 4.195:
        MCLmassLow[i] = 9.0
    else:
        MCLmassLow[i] = 6.91926521 + 0.49598529*ADx['redshift_BEAGLE'][i]
idx4 = (ADx['mass_BEAGLE_stellar'] + ADx['mag_AD'] > MCLmassLow)

idx5_z1 = (ADx['redshift_BEAGLE'] > z_lower) & (ADx['redshift_BEAGLE'] < z_upper)
idx5_z2 = (abs(ADx['redshift_BEAGLE']-ADx['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/visual_inspection/vis3/vis3_selection.csv', delimiter=",", skip_header=1)
idx5_z3 = np.full(len(ADx['id_AD']), False)
for i in range(len(ADx['id_AD'])):
    idx5_z3_temp = np.isin(ADx['id_AD'][i],vis[:,1][(vis[:,0]==ADx['field_AD'][i])&(vis[:,5]==1)])
    if idx5_z3_temp:
        idx5_z3[i] = True
idx5_z5 = (idx5_z3) | ((ADx['redshift_BEAGLE'] < 3.5) & (ADx['redshift_AD'] < 3.5) & (idx5_z1) & (idx5_z2))   

TmassHigh = 9.244 + (0.753*np.minimum(ADx['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(ADx['redshift_BEAGLE'], 4.0)**2)) # all galaxies
idx6 = (ADx['mass_BEAGLE_stellar'] < (TmassHigh)) 

idx7 = (ADx['min_chi2_BEAGLE']>0) & (ADx['min_chi2_BEAGLE']<13.28) # chi-squared
idx8 = (ADx['sfr_BEAGLE_instant']>-5.0) & (ADx['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
idx9 = (ADx['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103


# =============================================================================
# the plot (z vs z)
# =============================================================================

idx_contour = np.logical_and(idx2,idx3) #relflag + H
idx_contour = np.logical_and(idx_contour,idx7) #chi2
idx_KIRAC = np.logical_and(idx_contour,idx3_KIRAC_removed) # idx for no KIRAC objects
idx_subset = np.logical_and(idx_contour,idx3_KIRAC) # KIRAC cut
idx_subset = np.logical_and(idx_subset,idx1) # clusters at the mo
idx_subset = np.logical_and(idx_subset,idx4) # lower mass cut
idx_subset = np.logical_and(idx_subset,idx6) # upper mass cut
idx_subset = np.logical_and(idx_subset,idx5_z5) # redshift cuts + visual inspection
idx_subset = np.logical_and(idx_subset,idx8) # arbitrary sfr cut
idx_subset = np.logical_and(idx_subset,idx9) # GMM cut
print(sum(idx_subset))

ADx_subset = copy.deepcopy(ADx)
for key in ADx_subset.keys():
    ADx_subset[key] = ADx_subset[key][idx_subset]
    
    
###########
nbins_x = 30
nbins_y = 30
min_x = np.min(ADx['redshift_AD'][idx_contour])
max_x = np.max(ADx['redshift_AD'][idx_contour])
binsize_x = (max_x-min_x)/nbins_x
min_y = np.min(ADx['redshift_BEAGLE'][idx_contour])
max_y = np.max(ADx['redshift_BEAGLE'][idx_contour])
binsize_y = (min_y-max_y)/nbins_y
ADx1 = plt.hist2d(ADx['redshift_AD'][idx_contour], ADx['redshift_BEAGLE'][idx_contour], bins=[40,40], range=[[min_x,max_x],[min_y,max_y]])
###########

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
ax1.contourf((ADx1[1][1:]+ADx1[1][:-1])/2., \
                (ADx1[2][1:]+ADx1[2][:-1])/2., np.log10(np.transpose(ADx1[0])), \
                 cmap=cm.gist_yarg)
ax1.scatter(ADx['redshift_AD'][idx_KIRAC], ADx['redshift_BEAGLE'][idx_KIRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
ax1.scatter(ADx_subset['redshift_AD'], ADx_subset['redshift_BEAGLE'], s=3, alpha=1.0, color='red', label='Selected Subset')
#ax1.plot((0, 10),(0, 10), color='lime')

ax1.set_xlim(0.0, 10.0)
ax1.set_xlabel(r'ASTRODEEP redshift', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(0.0, 10.0)
ax1.set_ylabel(r'BEAGLE redshift', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

handles, labels = ax1.get_legend_handles_labels()
patch = mpatches.Patch(color='grey', label='Measurable Photo-z')
handles.append(patch) 
handles = [handles[2], handles[0], handles[1]]

ax1.legend(handles=handles, bbox_to_anchor=(0.95, 0.05), loc=4, frameon=True, fontsize=fontsize_legend)

if save:
    plt.savefig('002_redshift_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# redshift histograms
# =============================================================================

fig = plt.figure(figsize=(0.6*figuresize, 0.6*figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
#ax1.contourf((ADx1[1][1:]+ADx1[1][:-1])/2., \
#                (ADx1[2][1:]+ADx1[2][:-1])/2., np.log10(np.transpose(ADx1[0])), \
#                 cmap=cm.gist_yarg)
#ax1.scatter(ADx['redshift_AD'][idx_IRAC], ADx['redshift_BEAGLE'][idx_IRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
#ax1.scatter(ADx_subset['redshift_AD'], ADx_subset['redshift_BEAGLE'], s=3, alpha=1.0, color='red', label='Selected Subset')


ax1.hist(((ADx_subset['redshift_BEAGLE']-ADx_subset['redshift_AD'])/(1.0+ADx_subset['redshift_AD'])), alpha=0.3, bins=40,log=True, label='median z')
#ax1.hist(((ADx_subset['redshift_BEAGLE_mean']-ADx_subset['redshift_AD'])/(1.0+ADx_subset['redshift_AD'])), alpha=0.3, bins=40, log=True, label='mean z')

#ax1.hist(ADx_subset['redshift_AD'], alpha=0.3, bins=30)
#ax1.hist(ADx_subset['redshift_BEAGLE'], alpha=0.3, bins=30)

#ax1.set_xlim(0.0, 10.0)
ax1.set_xlabel(r'$(z_\mathrm{BEAGLE}-z_\mathrm{ASTRODEEP})/(1+z_\mathrm{ASTRODEEP})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.set_ylim(0.0, 10.0)
ax1.set_ylabel(r'Count', labelpad=10)
#ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)


ax1.legend(bbox_to_anchor=(0.95, 0.95), loc=1, frameon=True, fontsize=fontsize_legend)

if save:
    plt.savefig('010_redshift_histograms.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# lower mass limit vs redshift plot
# =============================================================================

mcl_mass_new_90 = mass_completeness_limits_0p90_new[1].data['mass']
mcl_z_new_90 = mass_completeness_limits_0p90_new[1].data['redshift']
fit_new_90 = np.polyfit(mcl_z_new_90, mcl_mass_new_90, 1)

'''
I think you could safely fit the first part with a straight line
 from z~2-4.5 or so, put a minimum at 8.1 and plateau at ~9.1 
 (wherever your fitted line crosses 9.1).
'''

fit_new_90_2_to_4p5 = np.polyfit(mcl_z_new_90[(mcl_z_new_90>2.0)&(mcl_z_new_90<4.5)], mcl_mass_new_90[(mcl_z_new_90>2.0)&(mcl_z_new_90<4.5)], 1)

mass_low_90 = 8.0
mass_high_90 = 9.0

z_low_90 = (mass_low_90 - fit_new_90_2_to_4p5[1]) / fit_new_90_2_to_4p5[0]
z_high_90 = (mass_high_90 - fit_new_90_2_to_4p5[1]) / fit_new_90_2_to_4p5[0]

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]

ax1.scatter(mcl_z_new_90, mcl_mass_new_90, color='k', marker='x', s=100, zorder=3, linewidth=3)
ax1.plot(np.linspace(z_low_90, z_high_90, 2), fit_new_90_2_to_4p5[1]+np.linspace(z_low_90, z_high_90, 2)*fit_new_90_2_to_4p5[0], color='grey', label='Lower Mass Limit', linewidth=3)
ax1.plot((0.0, z_low_90), (mass_low_90, mass_low_90), color='grey', linewidth=3)
ax1.plot((z_high_90, 8.0), (mass_high_90, mass_high_90), color='grey', linewidth=3)


ax1.scatter(ADx_subset['redshift_BEAGLE'], ADx_subset['mass_BEAGLE_stellar'] + ADx_subset['mag_AD'], s=3, alpha=1.0, color='blue', label='BEAGLE Fitted Mass', zorder=2)

ax1.scatter(ADx_subset['redshift_BEAGLE'], ADx_subset['mass_BEAGLE_stellar'], s=3, alpha=1.0, color='red', label='Magnification Corrected Mass', zorder=1)

ax1.set_xlim(0.0, 8.0)
ax1.set_xlabel(r'BEAGLE redshift', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(6.7, 11.3)
ax1.set_ylabel(string_mass, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

ax1.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=True, fontsize=fontsize_legend)

if save:
    plt.savefig('009_lower_mass_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()







#%%
# =============================================================================
# MS colour coded by redshift
# =============================================================================

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

#ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')

x_tmp = np.array((0, 20))
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+2)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 2*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+3)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 3*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+4)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 4*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+5)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 5*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+6)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 6*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')

ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

idx_sort = np.argsort(ADx_subset['redshift_BEAGLE'])
scatter = ax1.scatter(ADx_subset['mass_BEAGLE_stellar'][idx_sort], ADx_subset['sfr_BEAGLE_instant'][idx_sort], c=ADx_subset['redshift_BEAGLE'][idx_sort], vmin=1.25, vmax=6.5, alpha=1.0)

ax1.set_xlim(6.5, 10.5)
ax1.set_xlabel(string_mass, labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(string_sfr, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label('Redshift', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

if save:
    plt.savefig('003_MS_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# Delta MS colour coded by redshift
# =============================================================================

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

sfr_surface_real = ((medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b'])*ADx_subset['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+ADx_subset['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b']))

idx_sort = np.argsort(ADx_subset['redshift_BEAGLE'])

scatter = ax1.scatter(ADx_subset['mass_BEAGLE_stellar'][idx_sort], (ADx_subset['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=ADx_subset['redshift_BEAGLE'][idx_sort], vmin=1.25, vmax=6.5, alpha=1.0)

x_test = np.linspace(5, 12, 10)
ax1.plot(x_test, x_test*0, color='k', alpha=0.5)

ax1.set_xlim(5.5, 11.5)
ax1.set_xlabel(string_mass, labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(string_deltaMS, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)


cb.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label('Redshift', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

if save:
    plt.savefig('004_MS_collapsed_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# MS colour coded by MS probability
# =============================================================================

sfr_surface_real = ((medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b'])*ADx_subset['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+ADx_subset['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b']))

log_p_xi_eta_theta = norm.logpdf(ADx_subset['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(ADx_subset['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
p_bad = norm.pdf(ADx_subset['sfr_BEAGLE_instant'], scale=medians['outlier_sigma'], loc=medians['outlier_mean'])

z_bad = medians['pbad']*p_bad
z_good = (1-medians['pbad'])*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

#ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')

x_tmp = np.array((0, 20))
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+2)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 2*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+3)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 3*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+4)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 4*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+5)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 5*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(medians['alphaN_a']*(1+6)**medians['alphaN_b']) + normalisation - 9.0 + (medians['beta_a'] + 6*medians['beta_b'])*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')

ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

scatter = ax1.scatter(ADx_subset['mass_BEAGLE_stellar'][idx_sort], ADx_subset['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5)

ax1.set_xlim(6.5, 10.5)
ax1.set_xlabel(string_mass, labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(string_sfr, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(string_prob_ratio, rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

if save:
    plt.savefig('005_MS_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# Delta MS colour coded by MS probability
# =============================================================================

sfr_surface_real = ((medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b'])*ADx_subset['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+ADx_subset['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+ADx_subset['redshift_BEAGLE']*medians['beta_b']))

log_p_xi_eta_theta = norm.logpdf(ADx_subset['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(ADx_subset['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
p_bad = norm.pdf(ADx_subset['sfr_BEAGLE_instant'][idx_sort], scale=medians['outlier_sigma'], loc=medians['outlier_mean'])

z_bad = medians['pbad']*p_bad
z_good = (1-medians['pbad'])*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

scatter = ax1.scatter(ADx['mass_BEAGLE_stellar'][idx_sort], (ADx_subset['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5)

x_test = np.linspace(5, 12, 10)
ax1.plot(x_test, x_test*0, color='k', alpha=0.5)

ax1.set_xlim(5.5, 11.5)
ax1.set_xlabel(string_mass, labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(string_deltaMS, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(string_prob_ratio, rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

if save:
    plt.savefig('006_MS_collapsed_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# bias test histograms
# =============================================================================

names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
names_plot = ['alphaN a', 'alphaN b', 'beta a', 'beta b', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
trues = [4.0, -0.75, 0.5, 0.25, 0.3, 0.3, 0.0, 2.0] # outlier_mean is REALLY 0

dic = {}
for name in names:
    dic[name+'_16'] = []
    dic[name] = []
    dic[name+'_84'] = []
    
for i in range(len(bias_tests)): 
    nChains=4
    minIter=len(chain_bias[i])/nChains
    burn=int(0.75*minIter)
#        burn=0
    chain_arr = []
    for j in range(nChains):
        start = int(minIter*j+burn)
        finish = int(minIter*(j+1))
#        print(start, finish)
        chain_arr.append(chain_bias[i][start:finish])
    chain_combined = np.concatenate(chain_arr)   
    for name in names:
        dic[name+'_16'].append(float('{:.3f}'.format(np.percentile(chain_combined[name], 16))))
        dic[name].append(float('{:.3f}'.format(np.median(chain_combined[name]))))
        dic[name+'_84'].append(float('{:.3f}'.format(np.percentile(chain_combined[name], 84))))

#plot absolute around 0
x_values = np.linspace(1.5, 2., 10)
fig = plt.figure(figsize=(2*figuresize, 0.5*figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
for j, name in enumerate(names):
    ax1.scatter(x_values+j, np.array(dic[name]) - trues[j], color='k', marker='x')
    ax1.plot((x_values+j, x_values+j), (np.array(dic[name+'_16']) - trues[j], np.array(dic[name+'_84']) - trues[j]), color='k', zorder=0, alpha=0.5)   
ax1.plot((0,10),(0,0), color='k', alpha=0.5)

ax1.set_xlim(1.25, 9.25)
ax1.set_xticks(np.linspace(1.75, 8.75, len(names_plot)))
ax1.set_xticklabels(names_plot)
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-0.75, 0.75)
ax1.set_ylabel(string_bias_test, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

if save:
    plt.savefig('007_bias_test_histograms.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


sys.exit()




#%%
## =============================================================================
## Corner Plots
## =============================================================================

chain_corner = chain_MS_274
    
#names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
#names_plot = ['alphaN a', 'alphaN b', 'beta a', 'beta b', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
#data_corner = np.array([chain_corner['alphaN_a'],chain_corner['alphaN_b'],chain_corner['beta_a'],chain_corner['beta_b'],chain_corner['sig0'],chain_corner['pbad'],chain_corner['outlier_mean'],chain_corner['outlier_sigma']]).T
    
#names = ['alphaN_a', 'alphaN_b', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
data_corner = np.array([chain_corner['alphaN_a'],chain_corner['alphaN_b'],chain_corner['beta_a'],chain_corner['sig0'],chain_corner['pbad'],chain_corner['outlier_mean'],chain_corner['outlier_sigma']]).T

figure = corner.corner(data_corner, labels=names_plot,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": fontsize_axes})

if save:
    plt.savefig('008_corner_plots.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%% 
# =============================================================================
# heatplots from GMM, 1.25<z<2 with different lower mass cuts
# =============================================================================

num = 3
z_med_hp_low = 1.25
z_med_hp_high = 2.0
hp_lower_masses = [8.0, 8.5, 9.0]

z_med_hp = (z_med_hp_low+z_med_hp_high)/2.0
z_med_hp_gap = (z_med_hp_low+z_med_hp_high)/2.0 - z_med_hp_low
santini_idx = 1 # 1.3 to 2.0

for m in range(len(hp_lower_masses)):
    
    idx_rdm = np.arange(len(s23))[(s23['redshift_BEAGLE']>z_med_hp_low)&(s23['redshift_BEAGLE']<z_med_hp_high)&(s23['mass_BEAGLE_stellar']+s23['mag_AD']>hp_lower_masses[m])] 
    print(len(idx_rdm))
    
    
    
    plt.hist((s23['mass_BEAGLE_stellar']+s23['mag_AD'])[idx_rdm], histtype='step', bins=30)
    plt.hist(s23['mass_BEAGLE_stellar'][idx_rdm], histtype='step', bins=30)
    plt.show()
    
    
    
    plt.hist(s23['sfr_BEAGLE_instant'][idx_rdm], histtype='step', bins=30)
    plt.show()    
    
    
    
    x_hp = np.array([])
    y_hp = np.array([])
    z_hp = np.array([])
    
    n_hp = 300 # number of samples to take from GMM in total
    
    for i in idx_rdm:

        for G in range(3):
            
            mean = np.array([s23['x_GMM_3d'][i,G],s23['y_GMM_3d'][i,G],s23['z_GMM_3d'][i,G]])
            cov = np.array([[np.power(s23['xsig_GMM_3d'][i,G],2), s23['xycov_GMM_3d'][i,G], s23['xzcov_GMM_3d'][i,G]],[s23['xycov_GMM_3d'][i,G], np.power(s23['ysig_GMM_3d'][i,G],2), s23['yzcov_GMM_3d'][i,G]],[s23['xzcov_GMM_3d'][i,G], s23['yzcov_GMM_3d'][i,G], np.power(s23['zsig_GMM_3d'][i,G],2)]])
    
            xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s23['amp_GMM_3d'][i,G]))
    
            x_hp = np.concatenate((x_hp,xyz[:,0]))
            y_hp = np.concatenate((y_hp,xyz[:,1]))
            z_hp = np.concatenate((z_hp,xyz[:,2]))

    # only keep GMM samples within the redshift bin
    x_hp = x_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
    y_hp = y_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
    z_hp = z_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]

    fig = plt.figure(figsize=(0.7*figuresize, 0.7*figuresize))
    
    ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
#    ax1.set_title('{} - {}'.format(z_med_hp_low, z_med_hp_high))
    
    xlow = 5.5
    xhigh = 11.5
    ylow = -8.5
    yhigh = 3.5
    
    h = ax1.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
    ax1.plot((xlow,xhigh), (alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xlow,alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xhigh), color='w') # santini
    ax1.plot((9.7,9.7), (ylow, yhigh), color='w')

    ax1.scatter(s23['mass_BEAGLE_stellar'][idx_rdm], s23['sfr_BEAGLE_instant'][idx_rdm], marker='x', color='r')


    if m == 0:
        ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = ssfr_a_275,ssfr_b_275,beta_a_275,beta_b_275
    elif m == 1:
        ssfr_a_tmp, sfr_b_tmp,beta_a_tmp,beta_b_tmp = ssfr_a_276,ssfr_b_276,beta_a_276,beta_b_276
    elif m == 2:
        ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = ssfr_a_277,ssfr_b_277,beta_a_277,beta_b_277
    
#    redshift_tmp = np.linspace(z_med_hp_low, z_med_hp_high, num)
    redshift_tmp = z_med_hp
    beta_tmp = beta_a_tmp + redshift_tmp*beta_b_tmp
    alpha_tmp = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + normalisation - 9.0
    ssfr_tmp = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

    ax1.plot((xlow,xhigh), (alpha_tmp + beta_tmp*(xlow-normalisation), alpha_tmp + beta_tmp*(xhigh-normalisation)), color='r')

    ax1.set_xlim(xlow, xhigh)
    ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
    ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax1.xaxis.set_tick_params(labelsize=fontsize_axes)
    
    ax1.set_ylim(ylow, yhigh)
    ax1.set_ylabel(r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
    ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
    ax1.yaxis.set_tick_params(labelsize=fontsize_axes)
    
    cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
    cb = plt.colorbar(h[3], cax=cbaxes)
    cb.set_ticks([3, 5, 10, 20])
    cb.set_ticklabels([3, 5, 10, 20])
    #cbaxes.set_ylabel(r'TEST', rotation=270, labelpad=30)
    #cbaxes.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=1))
    cbaxes.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    cbaxes.yaxis.set_tick_params(which='minor', size=0)
    cbaxes.yaxis.set_tick_params(labelsize=fontsize_axes)
    
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax1.set_facecolor(rgba)
    
#    ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)
    
    if save:
        plt.savefig('011_heatplot_z1to2_lowermass{}.png'.format(str(hp_lower_masses[m]).replace('.','p')), dpi=300, transparent=False, bbox_inches='tight')
    plt.show()
    

#%%
# =============================================================================
# redshift histograms
# =============================================================================

idx_santini = ((ADx_subset['mass_BEAGLE_stellar']+ADx_subset['mag_AD']) > 9.3) & (ADx_subset['redshift_BEAGLE'] > 1.0) & (ADx_subset['redshift_BEAGLE'] < 2.0)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
#ax1.contourf((ADx1[1][1:]+ADx1[1][:-1])/2., \
#                (ADx1[2][1:]+ADx1[2][:-1])/2., np.log10(np.transpose(ADx1[0])), \
#                 cmap=cm.gist_yarg)
#ax1.scatter(ADx['redshift_AD'][idx_IRAC], ADx['redshift_BEAGLE'][idx_IRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
#ax1.scatter(ADx_subset['redshift_AD'], ADx_subset['redshift_BEAGLE'], s=3, alpha=1.0, color='red', label='Selected Subset')


ax1.hist((ADx_subset['mass_BEAGLE_stellar']-ADx_subset['mass_SANTINI'])[idx_santini], alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
ax1.hist((ADx_subset['sfr_BEAGLE_instant']-ADx_subset['sfr_SANTINI'])[idx_santini], alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])


#ax1.hist(((ADx_subset['mass_BEAGLE_stellar']-ADx_subset['mass_SANTINI'])/(1.0+ADx_subset['mass_SANTINI'])), alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
#ax1.hist(((ADx_subset['sfr_BEAGLE_instant']-ADx_subset['sfr_SANTINI'])/(1.0+ADx_subset['sfr_SANTINI'])), alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])

#ax1.hist(ADx_subset['redshift_AD'], alpha=0.3, bins=30)
#ax1.hist(ADx_subset['redshift_BEAGLE'], alpha=0.3, bins=30)

#ax1.set_xlim(-50, 50)
ax1.set_xlabel(r'$\mathrm{BEAGLE}-\mathrm{SANTINI}$', labelpad=10)
#ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.set_ylim(0.0, 10.0)
ax1.set_ylabel(r'Count', labelpad=10)
#ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)


ax1.legend(bbox_to_anchor=(0.95, 0.95), loc=1, frameon=True, fontsize=fontsize_legend)

if save:
    plt.savefig('013_beagle_santini_mass_sfr_hist.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%

# =============================================================================
# Parameter vs redshift, in redshift bins ORIGINAL
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = 0.3
xhigh = 5.7

# =============================================================================
# beta
# =============================================================================

ax1 = fig.add_axes([0, 2.5, 0.5, 0.5]) #[left, bottom, width, height]

# redshift bins
#ax1.plot(z267, beta_a_267 + z267*beta_b_267, color='blue') #1.0
ax1.plot(z268, beta_a_268 + z268*beta_b_268, color='blue') #2.0
ax1.plot(z269, beta_a_269 + z269*beta_b_269, color='blue') #3.0
ax1.plot(z270, beta_a_270 + z270*beta_b_270, color='blue') #4.0
ax1.plot(z273, beta_a_273 + z273*beta_b_273, color='blue') #1.0
#ax1.fill_between(z267, beta_16_arr_267, beta_84_arr_267, alpha=0.3, color='blue', zorder=0)
ax1.fill_between(z268, beta_16_arr_268, beta_84_arr_268, alpha=0.3, color='blue', zorder=0)
ax1.fill_between(z269, beta_16_arr_269, beta_84_arr_269, alpha=0.3, color='blue', zorder=0)
ax1.fill_between(z270, beta_16_arr_270, beta_84_arr_270, alpha=0.3, color='blue', zorder=0)
ax1.fill_between(z273, beta_16_arr_273, beta_84_arr_273, alpha=0.3, color='blue', zorder=0)

ax1.set_xlim(xlow, xhigh)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')

ax1.set_ylim(0.1, 1.3)
#ax1.set_ylim(0.3, 1.3)
ax1.set_ylabel(string_slope, labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# alpha
# =============================================================================

ax2 = fig.add_axes([0, 2.0, 0.5, 0.5]) #[left, bottom, width, height]

#print(np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0)

# redshift bins
#ax2.plot(z267, np.log10(ssfr_a_267*(1.0+z267)**ssfr_b_267) + normalisation - 9.0, color='blue') #1.0
ax2.plot(z268, np.log10(ssfr_a_268*(1.0+z268)**ssfr_b_268) + normalisation - 9.0, color='blue') #2.0
ax2.plot(z269, np.log10(ssfr_a_269*(1.0+z269)**ssfr_b_269) + normalisation - 9.0, color='blue') #3.0
ax2.plot(z270, np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0, color='blue', label='Redshift bins') #4.0
ax2.plot(z273, np.log10(ssfr_a_273*(1.0+z273)**ssfr_b_273) + normalisation - 9.0, color='blue') #1.0
#ax2.fill_between(z267, alpha_16_arr_267, alpha_84_arr_267, alpha=0.3, color='blue', zorder=0)
ax2.fill_between(z268, alpha_16_arr_268, alpha_84_arr_268, alpha=0.3, color='blue', zorder=0)
ax2.fill_between(z269, alpha_16_arr_269, alpha_84_arr_269, alpha=0.3, color='blue', zorder=0)
ax2.fill_between(z270, alpha_16_arr_270, alpha_84_arr_270, alpha=0.3, color='blue', zorder=0)
ax2.fill_between(z273, alpha_16_arr_273, alpha_84_arr_273, alpha=0.3, color='blue', zorder=0)

ax2.set_xlim(xlow, xhigh)
#ax2.set_xlabel('Redshift', labelpad=10)
ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax2.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
#ax2.xaxis.set_tick_params(labelsize=fontsize_axes)

ax2.set_ylim(0.1, 1.9)
ax2.set_ylabel(string_normalisation, labelpad=10)
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax2.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax2.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# scatter 
# =============================================================================

#sig0_267, sig0_16_267, sig0_84_267 = np.median(chain_MS_267['sig0']), np.percentile(chain_MS_267['sig0'], 16, axis=0), np.percentile(chain_MS_267['sig0'], 84, axis=0)
sig0_268, sig0_16_268, sig0_84_268 = np.median(chain_MS_268['sig0']), np.percentile(chain_MS_268['sig0'], 16, axis=0), np.percentile(chain_MS_268['sig0'], 84, axis=0)
sig0_269, sig0_16_269, sig0_84_269 = np.median(chain_MS_269['sig0']), np.percentile(chain_MS_269['sig0'], 16, axis=0), np.percentile(chain_MS_269['sig0'], 84, axis=0)
sig0_270, sig0_16_270, sig0_84_270 = np.median(chain_MS_270['sig0']), np.percentile(chain_MS_270['sig0'], 16, axis=0), np.percentile(chain_MS_270['sig0'], 84, axis=0)
sig0_273, sig0_16_273, sig0_84_273 = np.median(chain_MS_273['sig0']), np.percentile(chain_MS_273['sig0'], 16, axis=0), np.percentile(chain_MS_273['sig0'], 84, axis=0)

pbad_268, pbad_16_268, pbad_84_268 = np.median(chain_MS_268['pbad']), np.percentile(chain_MS_268['pbad'], 16, axis=0), np.percentile(chain_MS_268['pbad'], 84, axis=0)
pbad_269, pbad_16_269, pbad_84_269 = np.median(chain_MS_269['pbad']), np.percentile(chain_MS_269['pbad'], 16, axis=0), np.percentile(chain_MS_269['pbad'], 84, axis=0)
pbad_270, pbad_16_270, pbad_84_270 = np.median(chain_MS_270['pbad']), np.percentile(chain_MS_270['pbad'], 16, axis=0), np.percentile(chain_MS_270['pbad'], 84, axis=0)
pbad_273, pbad_16_273, pbad_84_273 = np.median(chain_MS_273['pbad']), np.percentile(chain_MS_273['pbad'], 16, axis=0), np.percentile(chain_MS_273['pbad'], 84, axis=0)

outlier_mean_268, outlier_mean_16_268, outlier_mean_84_268 = np.median(chain_MS_268['outlier_mean']), np.percentile(chain_MS_268['outlier_mean'], 16, axis=0), np.percentile(chain_MS_268['outlier_mean'], 84, axis=0)
outlier_mean_269, outlier_mean_16_269, outlier_mean_84_269 = np.median(chain_MS_269['outlier_mean']), np.percentile(chain_MS_269['outlier_mean'], 16, axis=0), np.percentile(chain_MS_269['outlier_mean'], 84, axis=0)
outlier_mean_270, outlier_mean_16_270, outlier_mean_84_270 = np.median(chain_MS_270['outlier_mean']), np.percentile(chain_MS_270['outlier_mean'], 16, axis=0), np.percentile(chain_MS_270['outlier_mean'], 84, axis=0)
outlier_mean_273, outlier_mean_16_273, outlier_mean_84_273 = np.median(chain_MS_273['outlier_mean']), np.percentile(chain_MS_273['outlier_mean'], 16, axis=0), np.percentile(chain_MS_273['outlier_mean'], 84, axis=0)

outlier_sigma_268, outlier_sigma_16_268, outlier_sigma_84_268 = np.median(chain_MS_268['outlier_sigma']), np.percentile(chain_MS_268['outlier_sigma'], 16, axis=0), np.percentile(chain_MS_268['outlier_sigma'], 84, axis=0)
outlier_sigma_269, outlier_sigma_16_269, outlier_sigma_84_269 = np.median(chain_MS_269['outlier_sigma']), np.percentile(chain_MS_269['outlier_sigma'], 16, axis=0), np.percentile(chain_MS_269['outlier_sigma'], 84, axis=0)
outlier_sigma_270, outlier_sigma_16_270, outlier_sigma_84_270 = np.median(chain_MS_270['outlier_sigma']), np.percentile(chain_MS_270['outlier_sigma'], 16, axis=0), np.percentile(chain_MS_270['outlier_sigma'], 84, axis=0)
outlier_sigma_273, outlier_sigma_16_273, outlier_sigma_84_273 = np.median(chain_MS_273['outlier_sigma']), np.percentile(chain_MS_273['outlier_sigma'], 16, axis=0), np.percentile(chain_MS_273['outlier_sigma'], 84, axis=0)
#print(np.median(chain_MS_267['pbad']), np.median(chain_MS_267['outlier_mean']), np.median(chain_MS_267['outlier_sigma']))
#print(np.median(chain_MS_268['pbad']), np.median(chain_MS_268['outlier_mean']), np.median(chain_MS_268['outlier_sigma']))
#print(np.median(chain_MS_269['pbad']), np.median(chain_MS_269['outlier_mean']), np.median(chain_MS_269['outlier_sigma']))
#print(np.median(chain_MS_270['pbad']), np.median(chain_MS_270['outlier_mean']), np.median(chain_MS_270['outlier_sigma']))
#
#print(len(chain_MS_267['pbad']))

ax3 = fig.add_axes([0, 1.5, 0.5, 0.5]) #[left, bottom, width, height]

# redshift bins
#ax3.plot(z267, np.full(len(z267), sig0_267), color='blue') #1.0
ax3.plot(z268, np.full(len(z268), sig0_268), color='blue') #2.0
ax3.plot(z269, np.full(len(z269), sig0_269), color='blue') #3.0
ax3.plot(z270, np.full(len(z270), sig0_270), color='blue') #4.0
ax3.plot(z273, np.full(len(z273), sig0_273), color='blue') #1.0
#ax3.fill_between(z267, np.full(len(z267), sig0_16_267), np.full(len(z267), sig0_84_267), alpha=0.3, color='blue', zorder=0)
ax3.fill_between(z268, np.full(len(z268), sig0_16_268), np.full(len(z268), sig0_84_268), alpha=0.3, color='blue', zorder=0)
ax3.fill_between(z269, np.full(len(z269), sig0_16_269), np.full(len(z269), sig0_84_269), alpha=0.3, color='blue', zorder=0)
ax3.fill_between(z270, np.full(len(z270), sig0_16_270), np.full(len(z270), sig0_84_270), alpha=0.3, color='blue', zorder=0)
ax3.fill_between(z273, np.full(len(z273), sig0_16_273), np.full(len(z273), sig0_84_273), alpha=0.3, color='blue', zorder=0)

ax3.set_xlim(xlow, xhigh)
#ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax3.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
#ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

ax3.set_ylim(0.0, 0.65)
#ax3.set_ylim(0.3, 1.3)
ax3.set_ylabel(string_scatter, labelpad=10)
ax3.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax3.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax3.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# pbad
# =============================================================================

ax4 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]

#print(np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0)

# redshift bins
#ax4.plot(z267, np.log10(ssfr_a_267*(1.0+z267)**ssfr_b_267) + normalisation - 9.0, color='blue') #1.0
ax4.plot(z268, np.full(len(z268), pbad_268), color='blue') #2.0
ax4.plot(z269, np.full(len(z269), pbad_269), color='blue') #3.0
ax4.plot(z270, np.full(len(z270), pbad_270), color='blue', label='Redshift bins') #4.0
ax4.plot(z273, np.full(len(z273), pbad_273), color='blue') #1.0
#ax4.fill_between(z267, alpha_16_arr_267, alpha_84_arr_267, alpha=0.3, color='blue', zorder=0)
ax4.fill_between(z268, np.full(len(z268), pbad_16_268), np.full(len(z268), pbad_84_268), alpha=0.3, color='blue', zorder=0)
ax4.fill_between(z269, np.full(len(z269), pbad_16_269), np.full(len(z269), pbad_84_269), alpha=0.3, color='blue', zorder=0)
ax4.fill_between(z270, np.full(len(z270), pbad_16_270), np.full(len(z270), pbad_84_270), alpha=0.3, color='blue', zorder=0)
ax4.fill_between(z273, np.full(len(z273), pbad_16_273), np.full(len(z273), pbad_84_273), alpha=0.3, color='blue', zorder=0)

ax4.set_xlim(xlow, xhigh)
#ax4.set_xlabel('Redshift', labelpad=10)
ax4.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax4.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax4.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax4.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
#ax4.xaxis.set_tick_params(labelsize=fontsize_axes)

ax4.set_ylim(0.0, 0.5)
ax4.set_ylabel(string_pbad, labelpad=10)
ax4.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax4.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax4.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax4.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax4.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax4.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# outlier mean
# =============================================================================

ax5 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]

#print(np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0)

# redshift bins
#ax5.plot(z267, np.log10(ssfr_a_267*(1.0+z267)**ssfr_b_267) + normalisation - 9.0, color='blue') #1.0
ax5.plot(z268, np.full(len(z268), outlier_mean_268), color='blue') #2.0
ax5.plot(z269, np.full(len(z269), outlier_mean_269), color='blue') #3.0
ax5.plot(z270, np.full(len(z270), outlier_mean_270), color='blue', label='Redshift bins') #4.0
ax5.plot(z273, np.full(len(z273), outlier_mean_273), color='blue') #1.0
#ax5.fill_between(z267, alpha_16_arr_267, alpha_84_arr_267, alpha=0.3, color='blue', zorder=0)
ax5.fill_between(z268, np.full(len(z268), outlier_mean_16_268), np.full(len(z268), outlier_mean_84_268), alpha=0.3, color='blue', zorder=0)
ax5.fill_between(z269, np.full(len(z269), outlier_mean_16_269), np.full(len(z269), outlier_mean_84_269), alpha=0.3, color='blue', zorder=0)
ax5.fill_between(z270, np.full(len(z270), outlier_mean_16_270), np.full(len(z270), outlier_mean_84_270), alpha=0.3, color='blue', zorder=0)
ax5.fill_between(z273, np.full(len(z273), outlier_mean_16_273), np.full(len(z273), outlier_mean_84_273), alpha=0.3, color='blue', zorder=0)

ax5.set_xlim(xlow, xhigh)
#ax5.set_xlabel('Redshift', labelpad=10)
ax5.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax5.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax5.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax5.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
#ax5.xaxis.set_tick_params(labelsize=fontsize_axes)

ax5.set_ylim(0.0, 1.5)
ax5.set_ylabel(string_outlier_mean, labelpad=10)
ax5.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax5.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax5.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax5.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax5.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax5.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# outlier sigma
# =============================================================================

ax6 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

#print(np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0)

# redshift bins
#ax6.plot(z267, np.log10(ssfr_a_267*(1.0+z267)**ssfr_b_267) + normalisation - 9.0, color='blue') #1.0
ax6.plot(z268, np.full(len(z268), outlier_sigma_268), color='blue') #2.0
ax6.plot(z269, np.full(len(z269), outlier_sigma_269), color='blue') #3.0
ax6.plot(z270, np.full(len(z270), outlier_sigma_270), color='blue', label='Redshift bins') #4.0
ax6.plot(z273, np.full(len(z273), outlier_sigma_273), color='blue') #1.0
#ax6.fill_between(z267, alpha_16_arr_267, alpha_84_arr_267, alpha=0.3, color='blue', zorder=0)
ax6.fill_between(z268, np.full(len(z268), outlier_sigma_16_268), np.full(len(z268), outlier_sigma_84_268), alpha=0.3, color='blue', zorder=0)
ax6.fill_between(z269, np.full(len(z269), outlier_sigma_16_269), np.full(len(z269), outlier_sigma_84_269), alpha=0.3, color='blue', zorder=0)
ax6.fill_between(z270, np.full(len(z270), outlier_sigma_16_270), np.full(len(z270), outlier_sigma_84_270), alpha=0.3, color='blue', zorder=0)
ax6.fill_between(z273, np.full(len(z273), outlier_sigma_16_273), np.full(len(z273), outlier_sigma_84_273), alpha=0.3, color='blue', zorder=0)

ax6.set_xlim(xlow, xhigh)
ax6.set_xlabel('Redshift', labelpad=10)
ax6.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax6.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax6.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax6.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax6.xaxis.set_tick_params(labelsize=fontsize_axes)

ax6.set_ylim(0.1, 1.9)
ax6.set_ylabel(string_outlier_sigma, labelpad=10)
ax6.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax6.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax6.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax6.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax6.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax6.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)










if save:
    plt.savefig('014_redshift_bins.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%




# =============================================================================
# mass dependent scatter heatplot 1.25<z<2.0
# =============================================================================



num = 3
z_med_hp_low = 1.25
z_med_hp_high = 2.0
hp_lower_masses = [8.0]

z_med_hp = (z_med_hp_low+z_med_hp_high)/2.0
z_med_hp_gap = (z_med_hp_low+z_med_hp_high)/2.0 - z_med_hp_low
santini_idx = 1 # 1.3 to 2.0

for m in range(len(hp_lower_masses)):
    
    idx_rdm = np.arange(len(s23))[(s23['redshift_BEAGLE']>z_med_hp_low)&(s23['redshift_BEAGLE']<z_med_hp_high)&(s23['mass_BEAGLE_stellar']+s23['mag_AD']>hp_lower_masses[m])] 
    print(len(idx_rdm))

    x_hp = np.array([])
    y_hp = np.array([])
    z_hp = np.array([])
    
    n_hp = 300 # number of samples to take from GMM in total
    
    for i in idx_rdm:

        for G in range(3):
            
            mean = np.array([s23['x_GMM_3d'][i,G],s23['y_GMM_3d'][i,G],s23['z_GMM_3d'][i,G]])
            cov = np.array([[np.power(s23['xsig_GMM_3d'][i,G],2), s23['xycov_GMM_3d'][i,G], s23['xzcov_GMM_3d'][i,G]],[s23['xycov_GMM_3d'][i,G], np.power(s23['ysig_GMM_3d'][i,G],2), s23['yzcov_GMM_3d'][i,G]],[s23['xzcov_GMM_3d'][i,G], s23['yzcov_GMM_3d'][i,G], np.power(s23['zsig_GMM_3d'][i,G],2)]])
    
            xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s23['amp_GMM_3d'][i,G]))
    
            x_hp = np.concatenate((x_hp,xyz[:,0]))
            y_hp = np.concatenate((y_hp,xyz[:,1]))
            z_hp = np.concatenate((z_hp,xyz[:,2]))

    # only keep GMM samples within the redshift bin
    x_hp = x_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
    y_hp = y_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
    z_hp = z_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]

    fig = plt.figure(figsize=(0.7*figuresize, 0.7*figuresize))
    
    ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
#    ax1.set_title('{} - {}'.format(z_med_hp_low, z_med_hp_high))
    
    xlow = 5.5
    xhigh = 11.5
    ylow = -3.5
    yhigh = 3.5
    
    ximin = 8.5
    ximax = 10.0
    

    
    h = ax1.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
#    ax1.plot((xlow,xhigh), (alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xlow,alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xhigh), color='w') # santini
    ax1.plot((9.7,9.7), (ylow, yhigh), color='w')

    if m == 0:
        ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = ssfr_a_278,ssfr_b_278,beta_a_278,beta_b_278
        sig0_tmp, k_tmp = np.median(chain_MS_278['sig0']), np.median(chain_MS_278['k'])

    
#    redshift_tmp = np.linspace(z_med_hp_low, z_med_hp_high, num)
    x_tmp = np.array((xlow, xhigh))
    redshift_tmp = z_med_hp
    beta_tmp = beta_a_tmp + redshift_tmp*beta_b_tmp
    alpha_tmp = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + normalisation - 9.0
    ssfr_tmp = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

    ax1.plot(x_tmp, alpha_tmp + beta_tmp*(x_tmp-normalisation), color='r')

    sig_tmp = sig0_tmp*(((1.0-k_tmp)*(x_tmp-ximax)/(ximax-ximin))+1.0)    
    ax1.plot(x_tmp, alpha_tmp + beta_tmp*(x_tmp-normalisation) + sig_tmp, color='r', linestyle='dashed')
    ax1.plot(x_tmp, alpha_tmp + beta_tmp*(x_tmp-normalisation) - sig_tmp, color='r', linestyle='dashed')

#8.5 and 10
    ax1.set_xlim(xlow, xhigh)
    ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
    ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
    ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
    ax1.xaxis.set_tick_params(labelsize=fontsize_axes)
    
    ax1.set_ylim(ylow, yhigh)
    ax1.set_ylabel(r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
    ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
    ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
    ax1.yaxis.set_tick_params(labelsize=fontsize_axes)
    
    cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
    cb = plt.colorbar(h[3], cax=cbaxes)
    cb.set_ticks([3, 5, 10, 20])
    cb.set_ticklabels([3, 5, 10, 20])
    #cbaxes.set_ylabel(r'TEST', rotation=270, labelpad=30)
    #cbaxes.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=1))
    cbaxes.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    cbaxes.yaxis.set_tick_params(which='minor', size=0)
    cbaxes.yaxis.set_tick_params(labelsize=fontsize_axes)
    
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax1.set_facecolor(rgba)
    
#    ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=fontsize_legend)
    
    if save:
        plt.savefig('015_heatplot_z1to2_scatter.png', dpi=300, transparent=False, bbox_inches='tight')
    plt.show()


#%%

# =============================================================================
# mass dependent scatter corner plot
# =============================================================================

chain_corner = chain_MS_278
    
#names_plot = ['ssfr a', 'beta a', 'sig0', 'k', 'pbad', 'outlier mean', 'outlier sigma']
#data_corner = np.array([chain_corner['alphaN_a'],chain_corner['beta_a'],chain_corner['sig0'],chain_corner['k'],chain_corner['pbad'],chain_corner['outlier_mean'],chain_corner['outlier_sigma']]).T

names_plot = ['ssfr a', 'beta a', 'sig0', 'k']
data_corner = np.array([chain_corner['alphaN_a'],chain_corner['beta_a'],chain_corner['sig0'],chain_corner['k']]).T

figure = corner.corner(data_corner, labels=names_plot,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": fontsize_axes})

if save:
    plt.savefig('016_heatplot_z1to2_corner.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%% 








