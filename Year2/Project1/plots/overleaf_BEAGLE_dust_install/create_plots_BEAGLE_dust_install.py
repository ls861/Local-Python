#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 08:41:18 2021

@author: lester
"""

import sys
sys.path.append(r'/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/site-packages')

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

from read_chain import read_chain


# for p in sys.path:
#     print(p)


# Collect all the font names available to matplotlib
#font_names = [f.name for f in fm.fontManager.ttflist]
#print(font_names)

# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
mpl.rc('image', cmap='jet')
cmap = mpl.cm.get_cmap('jet')
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['text.usetex'] = False
# plt.rcParams['text.usetex'] = True

figuresize = 7
fontsize_legend = 12
fontsize_axes = 20

# EAS POSTER
# fontsize_legend = 10
# fontsize_axes = 24

normalisation = 9.7

load = True
#load = False
# save = False
save = True





normalisation = 9.7

# =============================================================================
# useful axis strings
# =============================================================================

# string_slope = r'$\mathrm{MS\,Slope,}\,\beta$'
# string_normalisation = r'MS Normalisation, $\alpha$'
# string_scatter = r'$\mathrm{MS\,Scatter,}\,\sigma$'
# string_ssfr = r'$\log(\mathrm{sSFR}/\mathrm{yr}^{-1})$'
# string_mass = r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$'
# string_sfr = r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$'
# string_deltaMS = r'$\Delta_{MS}$'
# string_prob_ratio = r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$'
# string_bias_test = r'$\Delta \mathrm{Parameter}$'
# string_pbad = r'pbad'
# string_outlier_mean = r'outlier mean'
# string_outlier_sigma = r'outlier sigma'

# string_slope = r'$\mathrm{MS\,Slope,}\,\beta$'
# string_normalisation = r'MS Normalisation, $\alpha_{9.7}$'
# string_scatter = r'$\mathrm{MS\,Scatter,}\,\sigma$'
# string_pbad = r'pbad'
# string_outlier_mean = r'outlier mean'
# string_outlier_sigma = r'outlier sigma'

string_slope = r'$\beta$'
string_normalisation = r'$\alpha_\mathrm{9.7}$'
string_scatter = r'$\sigma$'
string_pbad = r'$p_\mathrm{OL}$'
string_outlier_mean = r'$\mu_\mathrm{OL}$'
string_outlier_sigma = r'$\sigma_\mathrm{OL}$'
string_ssfrNorm = r'$N$'
string_ssfrPower = r'$\gamma$'

# string_redshift = r'$z$'
string_redshift = r'$\mathrm{Redshift}$'

string_mass = r'$\log(M/ \, \mathrm{M_{\odot}})$'
string_sfr = r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$'
string_ssfr = r'$\log(\mathrm{sSFR} \, / \, \mathrm{yr}^{-1})$'

# string_mass = r'$\mathrm{M_{\star}}$'
# string_sfr = r'${\Psi}$'
# string_ssfr = r'${\Psi_\mathrm{s}}$'

string_deltaMS = r'$\Delta_{MS}$'
string_prob_ratio = r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$'
string_bias_test = r'$\Delta \mathrm{Parameter}$'


# =============================================================================
# opening data
# =============================================================================
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3

if load:

    
    folder = "/Users/lester/Documents/linmix_files/simple_chains/"

    # outlier model
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_26_005.p", 'rb') as f:
        chain_MS_33_26_c1 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_26_005.p", 'rb') as f:
        chain_MS_33_26_c2 = pickle.load(f, encoding='latin1') 
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_26_005.p", 'rb') as f:
        chain_MS_33_26_c3 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_26_005.p", 'rb') as f:
        chain_MS_33_26_c4 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_005.p", 'rb') as f:
        chain_MS_33_26_c5 = pickle.load(f, encoding='latin1')   
        
    # full run
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z1p25-6p0_4x20000_24_005.p", 'rb') as f:
        chain_MS_33_24_c = pickle.load(f, encoding='latin1')
    with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-6p0_4x20000_24_005.p", 'rb') as f:
        chain_MS_34_24_c = pickle.load(f, encoding='latin1')    
    
    
    # pbad = 0
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_28_005.p", 'rb') as f:
        chain_MS_33_28_c1 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_28_005.p", 'rb') as f:
        chain_MS_33_28_c2 = pickle.load(f, encoding='latin1') 
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_28_005.p", 'rb') as f:
        chain_MS_33_28_c3 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_28_005.p", 'rb') as f:
        chain_MS_33_28_c4 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_28_005.p", 'rb') as f:
        chain_MS_33_28_c5 = pickle.load(f, encoding='latin1')       
        
    # pbad = 0, sigma clipping
    with open(folder+"PROCESSED_lm_chain_scenario_35_clusters_z1p25-2p0_4x20000_28_006.p", 'rb') as f:
        chain_MS_35_28_c1 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_35_clusters_z2p0-3p0_4x20000_28_006.p", 'rb') as f:
        chain_MS_35_28_c2 = pickle.load(f, encoding='latin1') 
    with open(folder+"PROCESSED_lm_chain_scenario_35_clusters_z3p0-4p0_4x20000_28_006.p", 'rb') as f:
        chain_MS_35_28_c3 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_35_clusters_z4p0-5p0_4x20000_28_006.p", 'rb') as f:
        chain_MS_35_28_c4 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_35_clusters_z5p0-6p0_4x20000_28_006.p", 'rb') as f:
        chain_MS_35_28_c5 = pickle.load(f, encoding='latin1')        
    
    
    
    
    
    
    # folder = "/Users/lester/Documents/linmix_files/100000_sample_chains/"
        
    # # no clipping or outlier model
    # with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_19_001.p", 'rb') as f:
    #     chain_MS_29_19_c1 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_19_001.p", 'rb') as f:
    #     chain_MS_31_19_c2 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_19_001.p", 'rb') as f:
    #     chain_MS_31_19_c3 = pickle.load(f, encoding='latin1')          
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_19_001.p", 'rb') as f:
    #     chain_MS_31_19_c4 = pickle.load(f, encoding='latin1')          
    
    # # outlier model
    # with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_20_001.p", 'rb') as f:
    #     chain_MS_29_20_c1 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_20_001.p", 'rb') as f:
    #     chain_MS_29_20_c2 = pickle.load(f, encoding='latin1') 
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_20_001.p", 'rb') as f:
    #     chain_MS_31_20_c3 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p", 'rb') as f:
    #     chain_MS_31_20_c4 = pickle.load(f, encoding='latin1')  
        
    # # sigma clipping
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z1p25-2p0_4x50000_21_001.p", 'rb') as f:
    #     chain_MS_31_21_c1 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_21_001.p", 'rb') as f:
    #     chain_MS_31_21_c2 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_001.p", 'rb') as f:
    #     chain_MS_31_21_c3 = pickle.load(f, encoding='latin1')  
    # with open(folder+"PROCESSED_lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_21_001.p", 'rb') as f:
    #     chain_MS_31_21_c4 = pickle.load(f, encoding='latin1')

    # # z5-6 with outlier model
    # # with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_004.p", 'rb') as f:
    # #     chain_MS_33_26_c5 = pickle.load(f, encoding='latin1')    
    
    # # intrinsic scatter run
    # with open(folder+"PROCESSED_lm_chain_scenario_33_clusters_z1p25-2p0_4x5000_27_001.p", 'rb') as f:
    #     chain_MS_33_27_c1 = pickle.load(f, encoding='latin1')    
        
    # # full run
    # # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_24_002.p", 'rb') as f:
    # #     chain_MS_34_24_c = pickle.load(f, encoding='latin1')

    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-5p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c15 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-3p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c13 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z2p0-4p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c24 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z3p0-5p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c35 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-4p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c14 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z2p0-5p0_4x5000_24_003.p", 'rb') as f:
    #     chain_MS_34_24_c25 = pickle.load(f, encoding='latin1')
    # with open(folder+"PROCESSED_lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_22_001.p", 'rb') as f:
    #     chain_MS_34_22_c = pickle.load(f, encoding='latin1')



    # AD catalogue
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)

    # from /Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/replicating_santini_with_santini_input_3dGMM.py
    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p', 'rb') as f:
        ADx = pickle.load(f, encoding='latin1') 
    # with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_34_clusters_z1p25-6p0.p', 'rb') as f:
    #     ADx_subset = pickle.load(f, encoding='latin1') 
    #     #    print(ADx.keys())        





    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p90_new.fits'
    mass_completeness_limits_0p90_new = fits.open(fileName)
    #print(data_fits.info())
    #print(data_fits[1].header)

#     with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
#         data = pickle.load(f, encoding='latin1') 
# #    print(data.keys())

#     bias_tests = ['110', '111', '112', '113', '114', '115', '116', '117', '118', '119']
#     bias_tests = ['115', '117', '110', '119', '112', '111', '116', '113', '118', '114']
#     chain_bias = []
#     for test in bias_tests:
#         with open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_truncated_{}.p'.format(test), 'rb') as f:
#             chain_bias.append(pickle.load(f, encoding='latin1'))
            
    # fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.fits'
    # s23 = fits.open(fileName)
    # # print(s23.info())
    # # print(s23[1].header)
    # s23 = fits.open(fileName)[1].data # scenario 23

    folder = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/'
#    s29z1 = fits.open(fileName)
#    print(s29z1.info())
#    print(s29z1[1].header)

    fileNames = ['scenario_33_clusters_z1p25-2p0.fits',
                 'scenario_33_clusters_z2p0-3p0.fits',
                 'scenario_33_clusters_z3p0-4p0.fits',
                 'scenario_33_clusters_z4p0-5p0.fits',
                 'scenario_33_clusters_z5p0-6p0.fits',
                 'scenario_34_clusters_z1p25-6p0.fits']
    
    s33z1, s33z2, s33z3, s33z4, s33z5, s34 = fits.open(folder+fileNames[0])[1].data, fits.open(folder+fileNames[1])[1].data, fits.open(folder+fileNames[2])[1].data, fits.open(folder+fileNames[3])[1].data, fits.open(folder+fileNames[4])[1].data, fits.open(folder+fileNames[5])[1].data

    fileNames = ['scenario_31_clusters_z1p25-2p0.fits',
                 'scenario_31_clusters_z2p0-3p0.fits',
                 'scenario_31_clusters_z3p0-4p0.fits',
                 'scenario_31_clusters_z4p0-5p0.fits',
                 'scenario_31_clusters_z5p0-6p0.fits']
    
    s31z1, s31z2, s31z3, s31z4, s31z5 = fits.open(folder+fileNames[0])[1].data, fits.open(folder+fileNames[1])[1].data, fits.open(folder+fileNames[2])[1].data, fits.open(folder+fileNames[3])[1].data, fits.open(folder+fileNames[4])[1].data



folder = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/plots/site-packages/read_chain/'

if True:
    # # pbad = 0
    # rc_29_19_c1 = np.load(folder+'rc_29_19_c1.npy', allow_pickle=True)
    # rc_31_19_c2 = np.load(folder+'rc_31_19_c2.npy', allow_pickle=True)
    # rc_31_19_c3 = np.load(folder+'rc_31_19_c3.npy', allow_pickle=True)
    # rc_31_19_c4 = np.load(folder+'rc_31_19_c4.npy', allow_pickle=True)
    
    # # outlier model
    # rc_29_20_c1 = np.load(folder+'rc_29_20_c1.npy', allow_pickle=True)
    # rc_29_20_c2 = np.load(folder+'rc_29_20_c2.npy', allow_pickle=True)
    # rc_31_20_c3 = np.load(folder+'rc_31_20_c3.npy', allow_pickle=True)
    # rc_31_20_c4 = np.load(folder+'rc_31_20_c4.npy', allow_pickle=True)
    # # rc_33_26_c5 = np.load(folder+'rc_33_26_c5.npy', allow_pickle=True)       
    
    # # rc_33_26_c5 = read_chain(5.0, 6.0, chain_MS_33_26_c5, 'const_beta_ssfr_alpha_new')
    # # np.save(folder+'rc_33_26_c5.npy', rc_33_26_c5)

    
    # # sigma clipping
    # rc_31_21_c1 = np.load(folder+'rc_31_21_c1.npy', allow_pickle=True)
    # rc_31_21_c2 = np.load(folder+'rc_31_21_c2.npy', allow_pickle=True)
    # rc_31_21_c3 = np.load(folder+'rc_31_21_c3.npy', allow_pickle=True)
    # rc_31_21_c4 = np.load(folder+'rc_31_21_c4.npy', allow_pickle=True)
    
    # # full run w outlier model
    # # rc_34_24_c = np.load(folder+'rc_34_24_c.npy', allow_pickle=True)
    # rc_34_24_c15 = np.load(folder+'rc_34_24_c15.npy', allow_pickle=True)

    # # intrinsic scatter
    # rc_33_27_c1 = np.load(folder+'rc_33_27_c1.npy', allow_pickle=True)



    # # rc_34_24_c13 = read_chain(1.25, 3.0, chain_MS_34_24_c13, 'const_beta_ssfr_alpha_new')
    # # np.save(folder+'rc_34_24_c13.npy', rc_34_24_c13)
    # rc_34_24_c13 = np.load(folder+'rc_34_24_c13.npy', allow_pickle=True)

    # # rc_34_24_c24 = read_chain(2.0, 4.0, chain_MS_34_24_c24, 'const_beta_ssfr_alpha_new')
    # # np.save(folder+'rc_34_24_c24.npy', rc_34_24_c24)
    # rc_34_24_c24 = np.load(folder+'rc_34_24_c24.npy', allow_pickle=True)

    # rc_34_24_c35 = np.load(folder+'rc_34_24_c35.npy', allow_pickle=True)
    
    # # rc_34_24_c14 = read_chain(1.25, 4.0, chain_MS_34_24_c14, 'const_beta_ssfr_alpha_new')
    # # np.save(folder+'rc_34_24_c14.npy', rc_34_24_c14)
    # rc_34_24_c14 = np.load(folder+'rc_34_24_c14.npy', allow_pickle=True)

    # # rc_34_24_c25 = read_chain(2.0, 5.0, chain_MS_34_24_c25, 'const_beta_ssfr_alpha_new')
    # # np.save(folder+'rc_34_24_c25.npy', rc_34_24_c25)
    # rc_34_24_c25 = np.load(folder+'rc_34_24_c25.npy', allow_pickle=True)

    # rc_34_22_c = np.load(folder+'rc_34_22_c.npy', allow_pickle=True)



      
    
    # rc_33_26_c1 = read_chain(1.25, 2.0, chain_MS_33_26_c1, 'const_beta_linear_alpha')
    # rc_33_26_c2 = read_chain(2.0, 3.0, chain_MS_33_26_c2, 'const_beta_linear_alpha')
    # rc_33_26_c3 = read_chain(3.0, 4.0, chain_MS_33_26_c3, 'const_beta_linear_alpha')
    # rc_33_26_c4 = read_chain(4.0, 5.0, chain_MS_33_26_c4, 'const_beta_linear_alpha')
    # rc_33_26_c5 = read_chain(5.0, 6.0, chain_MS_33_26_c5, 'const_beta_linear_alpha')

    # rc_33_24_c = read_chain(1.25, 6.0, chain_MS_33_24_c, 'const_beta_ssfr_alpha_new')
    # rc_34_24_c = read_chain(1.25, 6.0, chain_MS_34_24_c, 'const_beta_ssfr_alpha_new')    
    
    
    # np.save(folder+'rc_33_26_c1.npy', rc_33_26_c1)
    # np.save(folder+'rc_33_26_c2.npy', rc_33_26_c2)
    # np.save(folder+'rc_33_26_c3.npy', rc_33_26_c3)
    # np.save(folder+'rc_33_26_c4.npy', rc_33_26_c4)
    # np.save(folder+'rc_33_26_c5.npy', rc_33_26_c5)
    
    # np.save(folder+'rc_33_24_c.npy', rc_33_24_c)
    # np.save(folder+'rc_34_24_c.npy', rc_34_24_c)


    # redshift bins
    rc_33_26_c1 = np.load(folder+'rc_33_26_c1.npy', allow_pickle=True) 
    rc_33_26_c2 = np.load(folder+'rc_33_26_c2.npy', allow_pickle=True) 
    rc_33_26_c3 = np.load(folder+'rc_33_26_c3.npy', allow_pickle=True) 
    rc_33_26_c4 = np.load(folder+'rc_33_26_c4.npy', allow_pickle=True) 
    rc_33_26_c5 = np.load(folder+'rc_33_26_c5.npy', allow_pickle=True) 
    
    # full run with and without Pick a Peak
    rc_33_24_c = np.load(folder+'rc_33_24_c.npy', allow_pickle=True) 
    rc_34_24_c = np.load(folder+'rc_34_24_c.npy', allow_pickle=True) 





    # rc_33_28_c1 = read_chain(1.25, 2.0, chain_MS_33_28_c1, 'const_beta_linear_alpha')
    # rc_33_28_c2 = read_chain(2.0, 3.0, chain_MS_33_28_c2, 'const_beta_linear_alpha')
    # rc_33_28_c3 = read_chain(3.0, 4.0, chain_MS_33_28_c3, 'const_beta_linear_alpha')
    # rc_33_28_c4 = read_chain(4.0, 5.0, chain_MS_33_28_c4, 'const_beta_linear_alpha')
    # rc_33_28_c5 = read_chain(5.0, 6.0, chain_MS_33_28_c5, 'const_beta_linear_alpha')
    
    # rc_35_28_c1 = read_chain(1.25, 2.0, chain_MS_35_28_c1, 'const_beta_linear_alpha')
    # rc_35_28_c2 = read_chain(2.0, 3.0, chain_MS_35_28_c2, 'const_beta_linear_alpha')
    # rc_35_28_c3 = read_chain(3.0, 4.0, chain_MS_35_28_c3, 'const_beta_linear_alpha')
    # rc_35_28_c4 = read_chain(4.0, 5.0, chain_MS_35_28_c4, 'const_beta_linear_alpha')
    # rc_35_28_c5 = read_chain(5.0, 6.0, chain_MS_35_28_c5, 'const_beta_linear_alpha')
    
    # np.save(folder+'rc_33_28_c1.npy', rc_33_28_c1)
    # np.save(folder+'rc_33_28_c2.npy', rc_33_28_c2)
    # np.save(folder+'rc_33_28_c3.npy', rc_33_28_c3)
    # np.save(folder+'rc_33_28_c4.npy', rc_33_28_c4)
    # np.save(folder+'rc_33_28_c5.npy', rc_33_28_c5)
    
    # np.save(folder+'rc_35_28_c1.npy', rc_35_28_c1)
    # np.save(folder+'rc_35_28_c2.npy', rc_35_28_c2)
    # np.save(folder+'rc_35_28_c3.npy', rc_35_28_c3)
    # np.save(folder+'rc_35_28_c4.npy', rc_35_28_c4)
    # np.save(folder+'rc_35_28_c5.npy', rc_35_28_c5)
    
    # pbad == 0
    rc_33_28_c1 = np.load(folder+'rc_33_28_c1.npy', allow_pickle=True) 
    rc_33_28_c2 = np.load(folder+'rc_33_28_c2.npy', allow_pickle=True) 
    rc_33_28_c3 = np.load(folder+'rc_33_28_c3.npy', allow_pickle=True) 
    rc_33_28_c4 = np.load(folder+'rc_33_28_c4.npy', allow_pickle=True) 
    rc_33_28_c5 = np.load(folder+'rc_33_28_c5.npy', allow_pickle=True) 
    
    # sigma clipping, pbad == 0
    rc_35_28_c1 = np.load(folder+'rc_35_28_c1.npy', allow_pickle=True) 
    rc_35_28_c2 = np.load(folder+'rc_35_28_c2.npy', allow_pickle=True) 
    rc_35_28_c3 = np.load(folder+'rc_35_28_c3.npy', allow_pickle=True) 
    rc_35_28_c4 = np.load(folder+'rc_35_28_c4.npy', allow_pickle=True) 
    rc_35_28_c5 = np.load(folder+'rc_35_28_c5.npy', allow_pickle=True) 






else:
    # outlier model
    rc_29_19_c1 = read_chain(1.25, 2.0, chain_MS_29_19_c1, 'const_beta_linear_alpha')
    rc_31_19_c2 = read_chain(2.0, 3.0, chain_MS_31_19_c2, 'const_beta_linear_alpha')
    rc_31_19_c3 = read_chain(3.0, 4.0, chain_MS_31_19_c3, 'const_beta_linear_alpha')
    rc_31_19_c4 = read_chain(4.0, 5.0, chain_MS_31_19_c4, 'const_beta_linear_alpha')
    rc_33_26_c5 = read_chain(5.0, 6.0, chain_MS_33_26_c5, 'const_beta_linear_alpha')
    
    # pbad = 0
    rc_29_20_c1 = read_chain(1.25, 2.0, chain_MS_29_20_c1, 'const_beta_linear_alpha')
    rc_29_20_c2 = read_chain(2.0, 3.0, chain_MS_29_20_c2, 'const_beta_linear_alpha')
    rc_31_20_c3 = read_chain(3.0, 4.0, chain_MS_31_20_c3, 'const_beta_linear_alpha')
    rc_31_20_c4 = read_chain(4.0, 5.0, chain_MS_31_20_c4, 'const_beta_linear_alpha')
    
    # sigma clipping
    rc_31_21_c1 = read_chain(1.25, 2.0, chain_MS_31_21_c1, 'const_beta_linear_alpha')
    rc_31_21_c2 = read_chain(2.0, 3.0, chain_MS_31_21_c2, 'const_beta_linear_alpha')
    rc_31_21_c3 = read_chain(3.0, 4.0, chain_MS_31_21_c3, 'const_beta_linear_alpha')
    rc_31_21_c4 = read_chain(4.0, 5.0, chain_MS_31_21_c4, 'const_beta_linear_alpha')   

    # full run w outlier model
    rc_34_24_c = read_chain(1.25, 6.0, chain_MS_34_24_c, 'const_beta_ssfr_alpha_new')
    rc_34_24_c15 = read_chain(1.25, 5.0, chain_MS_34_24_c15, 'const_beta_ssfr_alpha_new')
    rc_34_24_c35 = read_chain(3.0, 5.0, chain_MS_34_24_c35, 'const_beta_ssfr_alpha_new')
    
    rc_34_22_c = read_chain(1.25, 6.0, chain_MS_34_22_c, 'const_beta_ssfr_alpha')
    np.save(folder+'rc_34_22_c.npy', rc_34_22_c)

    # intrinsic scatter
    rc_33_27_c1 = read_chain(1.25, 2.0, chain_MS_33_27_c1, 'const_beta_linear_alpha')


    # SAVE
    np.save(folder+'rc_29_19_c1.npy', rc_29_19_c1)
    np.save(folder+'rc_31_19_c2.npy', rc_31_19_c2)
    np.save(folder+'rc_31_19_c3.npy', rc_31_19_c3)
    np.save(folder+'rc_31_19_c4.npy', rc_31_19_c4)
    np.save(folder+'rc_33_26_c5.npy', rc_33_26_c5)
    
    np.save(folder+'rc_29_20_c1.npy', rc_29_20_c1)
    np.save(folder+'rc_29_20_c2.npy', rc_29_20_c2)
    np.save(folder+'rc_31_20_c3.npy', rc_31_20_c3)
    np.save(folder+'rc_31_20_c4.npy', rc_31_20_c4)
    
    np.save(folder+'rc_31_21_c1.npy', rc_31_21_c1)
    np.save(folder+'rc_31_21_c2.npy', rc_31_21_c2)
    np.save(folder+'rc_31_21_c3.npy', rc_31_21_c3)
    np.save(folder+'rc_31_21_c4.npy', rc_31_21_c4)    

    np.save(folder+'rc_34_24_c.npy', rc_34_24_c)
    np.save(folder+'rc_34_24_c15.npy', rc_34_24_c15)
    np.save(folder+'rc_34_24_c35.npy', rc_34_24_c35)

    np.save(folder+'rc_33_27_c1.npy', rc_33_27_c1)   
    
    
# =============================================================================
# get medians for MS plots - eventually full run
# =============================================================================

def get_medians(chain_MS):
#    names = chain_original.dtype.names
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
    dic = {}
    for name in names:
        dic[name] = np.median(chain_MS[name])
    return dic

medians_34_24_c = get_medians(chain_MS_34_24_c)
# medians_33_27_c1 = get_medians(chain_MS_33_27_c1)



#%%

from redshift_bins_1 import redshift_bins_1
from redshift_bins_2 import redshift_bins_2
from redshift_bins_3 import redshift_bins_3
from redshift_bins_4 import redshift_bins_4
from redshift_bins_5 import redshift_bins_5
from full_run import full_run
from full_run_2 import full_run_2
from full_run_3 import full_run_3
from corner_full_run import corner_full_run
from redshift_vs_redshift import redshift_vs_redshift
from redshift_histogram import redshift_histogram
from lower_mass_limit import lower_mass_limit
from MS_redshift import MS_redshift
from deltaMS_redshift import deltaMS_redshift
from MS_probability import MS_probability
from deltaMS_probability import deltaMS_probability
from corner_scatter import corner_scatter
from heatplot_scatter import heatplot_scatter
from redshift_bins_MS import redshift_bins_MS

# redshift_bins_1(rc_29_19_c1, rc_31_19_c2, rc_31_19_c3, rc_31_19_c4, rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_31_21_c1, rc_31_21_c2, rc_31_21_c3, rc_31_21_c4, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)


# OLD RUNS
# redshift_bins_2(rc_29_19_c1, rc_31_19_c2, rc_31_19_c3, rc_31_19_c4, rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_33_26_c5, rc_31_21_c1, rc_31_21_c2, rc_31_21_c3, rc_31_21_c4, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

# redshift_bins_3(rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_33_26_c5, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

# redshift_bins_2(rc_29_19_c1, rc_31_19_c2, rc_31_19_c3, rc_31_19_c4, rc_33_26_c1, rc_33_26_c2, rc_33_26_c3, rc_33_26_c4, rc_33_26_c5, rc_31_21_c1, rc_31_21_c2, rc_31_21_c3, rc_31_21_c4, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

# NEW RUNS





# redshift_bins_4(rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, normalisation, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_ssfr, string_redshift, save)

# redshift_bins_5(rc_29_19_c1, rc_31_19_c2, rc_31_19_c3, rc_31_19_c4, rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_31_21_c1, rc_31_21_c2, rc_31_21_c3, rc_31_21_c4, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

# full_run(rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_34_24_c, normalisation, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_ssfr, save)

# full_run(rc_34_24_c15, rc_34_24_c15, rc_34_22_c, rc_34_22_c, rc_34_24_c, normalisation, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_ssfr, save)

# full_run(rc_34_24_c13, rc_34_24_c24, rc_34_24_c35, rc_34_24_c13, rc_34_24_c, normalisation, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_ssfr, save)

# full_run(rc_34_24_c14, rc_34_24_c14, rc_34_24_c25, rc_34_24_c25, rc_34_24_c, normalisation, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_ssfr, save)



# EAS POSTER
# fontsize_legend = 10
# fontsize_axes = 24

# entire code could do with fixing
# full_run_2(chain_MS_29_20_c1, chain_MS_29_20_c2, chain_MS_31_20_c3, chain_MS_31_20_c4, chain_MS_33_26_c5, chain_MS_33_27_c1, s33z1, s33z2, s33z3, s33z4, s33z5, rc_29_20_c1, rc_29_20_c2, rc_31_20_c3, rc_31_20_c4, rc_33_26_c5, rc_33_27_c1, rc_34_24_c, normalisation, figuresize, 10, 24, string_slope, string_normalisation, string_scatter, string_ssfr, save)


# save=True
# save=False




# redshift_histogram(s34, figuresize, fontsize_legend, fontsize_axes, save)




# corner_scatter(chain_MS_33_27_c1, fontsize_axes, save)
# heatplot_scatter(s33z1, medians_33_27_c1, normalisation, figuresize, fontsize_axes, string_mass, string_sfr, save)



# NO PICK A PEAK
# full_run_3(chain_MS_33_26_c1, chain_MS_33_26_c2, chain_MS_33_26_c3, chain_MS_33_26_c4, chain_MS_33_26_c5, s33z1, s33z2, s33z3, s33z4, s33z5, rc_33_26_c1, rc_33_26_c2, rc_33_26_c3, rc_33_26_c4, rc_33_26_c5, rc_33_24_c, normalisation, figuresize, 10, fontsize_axes, string_slope, string_normalisation, string_scatter, string_ssfr, string_redshift, save)

#%%


# fig 1
import create_filter_plot

# fig 3
redshift_vs_redshift(AD, ADx, s34, figuresize, 14, fontsize_axes, save)

# fig 4
lower_mass_limit(mass_completeness_limits_0p90_new, s34, figuresize, 14, fontsize_axes, string_mass, save)

# fig 5
redshift_bins_MS(chain_MS_33_26_c1, 1, 1.25, 2.0, normalisation, figuresize, fontsize_legend, fontsize_axes, string_mass, string_sfr, string_prob_ratio, False, save)

redshift_bins_MS(chain_MS_33_26_c2, 2, 2.0, 3.0, normalisation, figuresize, fontsize_legend, fontsize_axes, string_mass, string_sfr, string_prob_ratio, False, save)

redshift_bins_MS(chain_MS_33_26_c3, 3, 3.0, 4.0, normalisation, figuresize, fontsize_legend, fontsize_axes, string_mass, string_sfr, string_prob_ratio, False, save)

redshift_bins_MS(chain_MS_33_26_c4, 4, 4.0, 5.0, normalisation, figuresize, fontsize_legend, fontsize_axes, string_mass, string_sfr, string_prob_ratio, False, save)

redshift_bins_MS(chain_MS_33_26_c5, 5, 5.0, 6.0, normalisation, figuresize, 14, fontsize_axes, string_mass, string_sfr, string_prob_ratio, True, save)

# fig 6
redshift_bins_2(rc_33_28_c1, rc_33_28_c2, rc_33_28_c3, rc_33_28_c4, rc_33_28_c5, rc_33_26_c1, rc_33_26_c2, rc_33_26_c3, rc_33_26_c4, rc_33_26_c5, rc_35_28_c1, rc_35_28_c2, rc_35_28_c3, rc_35_28_c4, rc_35_28_c5, figuresize, 14, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

redshift_bins_3(rc_33_26_c1, rc_33_26_c2, rc_33_26_c3, rc_33_26_c4, rc_33_26_c5, figuresize, fontsize_legend, fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_redshift, save)

# fig 7
MS_redshift(s34, medians_34_24_c, normalisation, figuresize, fontsize_axes, string_mass, string_sfr, save)

deltaMS_redshift(s34, medians_34_24_c, normalisation, figuresize, fontsize_axes, string_mass, string_deltaMS, save)

MS_probability(s34, medians_34_24_c, normalisation, figuresize, fontsize_axes, string_mass, string_sfr, string_prob_ratio, save)

deltaMS_probability(s34, medians_34_24_c, normalisation, figuresize, fontsize_axes, string_mass, string_deltaMS, string_prob_ratio, save)

# fig 8 and 10
full_run_3(chain_MS_33_26_c1, chain_MS_33_26_c2, chain_MS_33_26_c3, chain_MS_33_26_c4, chain_MS_33_26_c5, s33z1, s33z2, s33z3, s33z4, s33z5, rc_33_26_c1, rc_33_26_c2, rc_33_26_c3, rc_33_26_c4, rc_33_26_c5, rc_34_24_c, normalisation, figuresize, 10, fontsize_axes, string_slope, string_normalisation, string_scatter, string_ssfr, string_redshift, save)

# fig 9
corner_full_run(chain_MS_34_24_c, '24', fontsize_axes, string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma, string_ssfrNorm, string_ssfrPower, save)

# fig 11
import create_prior_heatplot

# fig 12 and 13
import create_scatter_plots

# fig 14
import create_panels

# fig 15
import create_emission_lines_plot












#%%
# =============================================================================
# results tables in paper
# =============================================================================

# print(rc_34_24_c.item().keys())
# dict_keys(['z', 'beta_16_arr', 'beta_50_arr', 'beta_84_arr', 'ssfr_16_arr', 'ssfr_50_arr', 'ssfr_84_arr', 'alpha_16_arr', 'alpha_50_arr', 'alpha_84_arr', 'sig0_16_arr', 'sig0_50_arr', 'sig0_84_arr', 'k_16_arr', 'k_50_arr', 'k_84_arr', 'pbad_16_arr', 'pbad_50_arr', 'pbad_84_arr', 'outlier_mean_16_arr', 'outlier_mean_50_arr', 'outlier_mean_84_arr', 'outlier_sigma_16_arr', 'outlier_sigma_50_arr', 'outlier_sigma_84_arr'])


str1 = ['alpha','beta','sig0','pbad','outlier_mean','outlier_sigma']
str2 = ['\intercept','\slope','\scatter','\OLprob','\OLmean','\OLscatter']

for i in range(len(str1)):
    print('{} & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$  & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$  & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$  & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$  & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ \\\\'.format(str2[i],
                                                                                                                                                  
min(rc_33_26_c1.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c1.item()['{}_84_arr'.format(str1[i])]) - min(rc_33_26_c1.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c1.item()['{}_50_arr'.format(str1[i])]) - min(rc_33_26_c1.item()['{}_16_arr'.format(str1[i])]),                     

min(rc_33_26_c2.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c2.item()['{}_84_arr'.format(str1[i])]) - min(rc_33_26_c2.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c2.item()['{}_50_arr'.format(str1[i])]) - min(rc_33_26_c2.item()['{}_16_arr'.format(str1[i])]),    

min(rc_33_26_c3.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c3.item()['{}_84_arr'.format(str1[i])]) - min(rc_33_26_c3.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c3.item()['{}_50_arr'.format(str1[i])]) - min(rc_33_26_c3.item()['{}_16_arr'.format(str1[i])]),    

min(rc_33_26_c4.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c4.item()['{}_84_arr'.format(str1[i])]) - min(rc_33_26_c4.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c4.item()['{}_50_arr'.format(str1[i])]) - min(rc_33_26_c4.item()['{}_16_arr'.format(str1[i])]),    

min(rc_33_26_c5.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c5.item()['{}_84_arr'.format(str1[i])]) - min(rc_33_26_c5.item()['{}_50_arr'.format(str1[i])]),
min(rc_33_26_c5.item()['{}_50_arr'.format(str1[i])]) - min(rc_33_26_c5.item()['{}_16_arr'.format(str1[i])])

))



# normalisation for full run is redshift dependent
# ssfr_arr = ssfr_a_arr + ssfr_b_arr*np.log10(1.0+redshift_arr) - 9.0
print('')
# print(chain_MS_34_24_c.keys())
str3 = ['alphaN_a','alphaN_b','beta_a','sig0','pbad','outlier_mean','outlier_sigma']
str4 = ['\ssfrNorm','\ssfrPower','\slope','\scatter','\OLprob','\OLmean','\OLscatter']

for i in range(len(str3)):
    
    if str3[i]=='alphaN_a':
        print('{} & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ \\\\'.format(str4[i],
        np.median(10**chain_MS_34_24_c[str3[i]]),
        np.percentile(10**chain_MS_34_24_c[str3[i]], 84, axis=0) - np.median(10**chain_MS_34_24_c[str3[i]]),
        np.median(10**chain_MS_34_24_c[str3[i]]) - np.percentile(10**chain_MS_34_24_c[str3[i]], 16, axis=0)
        
        ))
        
        # print('{} & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$'.format(str4[i],
        # 10**np.median(chain_MS_34_24_c[str3[i]]),
        # 10**np.percentile(chain_MS_34_24_c[str3[i]], 84, axis=0) - 10**np.median(chain_MS_34_24_c[str3[i]]),
        # 10**np.median(chain_MS_34_24_c[str3[i]]) - 10**np.percentile(chain_MS_34_24_c[str3[i]], 16, axis=0)
        
        # ))

    else:
        print('{} & ${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$ \\\\'.format(str4[i],
        np.median(chain_MS_34_24_c[str3[i]]),
        np.percentile(chain_MS_34_24_c[str3[i]], 84, axis=0) - np.median(chain_MS_34_24_c[str3[i]]),
        np.median(chain_MS_34_24_c[str3[i]]) - np.percentile(chain_MS_34_24_c[str3[i]], 16, axis=0)
        
        ))  



'''
print(rc_34_24_c.item().keys())

N = 10**np.median(chain_MS_34_24_c['alphaN_a'])
g = np.median(chain_MS_34_24_c['alphaN_b'])

print(N, g)

z = rc_34_24_c.item()['z'][0]
ssfr = np.log10(N*((1+z)**g)) - 9.0
print(ssfr, rc_34_24_c.item()['ssfr_50_arr'][0])
# the difference is one is median of ssfr a, and median ssfr b, other is median of array of precalculated ssfrs (which is the correct way):
# rc file: np.median(ssfr_a_arr + ssfr_b_arr*np.log10(1.0+z) - 9.0)
# medians too early: np.median(ssfr_a_arr) + np.median(ssfr_b_arr)*np.log10(1.0+z) - 9.0
'''


#%%

# normalisation = 9.7

# z_arr = np.linspace(1.25, 6.0, 1000)

# chain_MS = chain_MS_34_24_c


# redshift_arr = np.repeat(np.array([z_arr]).T, len(chain_MS['beta_a']), axis=1).T

# beta_a_arr = np.array([chain_MS['beta_a']]).T
# beta_b_arr = np.array([chain_MS['beta_b']]).T
# #    print(np.shape(redshift_arr), np.shape(beta_b_arr))
# beta_arr = beta_a_arr + redshift_arr*beta_b_arr
# beta_16_arr = np.percentile(beta_arr, 16, axis=0)
# beta_50_arr = np.median(beta_arr, axis=0)
# beta_84_arr = np.percentile(beta_arr, 84, axis=0)

# # alphaN_a + alphaN_b*np.log10(1.0+z) + self.alphaNorm - 9.0 # new alphaN_a = old np.log10(alphaN_a)

# ssfr_a_arr = np.array([chain_MS['alphaN_a']]).T
# ssfr_b_arr = np.array([chain_MS['alphaN_b']]).T
# ssfr_arr = ssfr_a_arr + ssfr_b_arr*np.log10(1.0+redshift_arr) - 9.0
# ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
# ssfr_50_arr = np.median(ssfr_arr, axis=0)
# ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)

# alpha_arr = ssfr_arr + normalisation
# alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
# alpha_50_arr = np.median(alpha_arr, axis=0)
# alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)

# print(ssfr_50_arr)

# print(10**np.median(ssfr_a_arr), N)
# print(np.median(ssfr_b_arr), g)


# ssfr1 = np.log10(N*((1+z_arr)**g)) - 9.0

# ssfr1 = np.median(ssfr_a_arr + ssfr_b_arr*np.log10(1.0+z) - 9.0)
# ssfr2 = np.median(ssfr_a_arr) + np.median(ssfr_b_arr)*np.log10(1.0+z) - 9.0

# print(ssfr1, ssfr2)




#%%





#%%
# =============================================================================
# finding out how many from s34 have k c1 or c2
# =============================================================================
# # from final subset
# k = s34['Ks_BEAGLE_input']
# c1 = s34['CH1_BEAGLE_input']
# c2 = s34['CH2_BEAGLE_input']

# # idx = (s34['redshift_BEAGLE']>4.0) & (s34['redshift_BEAGLE']<4.5)
# # idx = s34['redshift_BEAGLE']>4.0
# idx = s34['redshift_BEAGLE']>4.5


# print(len(k[idx]))
# print('')
# print(sum(k[idx]>-60))
# print(sum(c1[idx]>-60))
# print(sum(c2[idx]>-60))
# print(sum((k[idx]>-60)|(c1[idx]>-60)|(c2[idx]>-60)))
# print('')
# print(sum((k[idx]>-60)&(c1[idx]>-60)))
# print(sum((k[idx]>-60)&(c2[idx]>-60)))
# print(sum((c1[idx]>-60)&(c2[idx]>-60)))
# print(sum(((k[idx]>-60)&(c1[idx]>-60))|((k[idx]>-60)&(c2[idx]>-60))|((c1[idx]>-60)&(c2[idx]>-60))))
# print('')
# print(sum(((k[idx]>-60)&(c1[idx]>-60))&(c2[idx]>-60)))



#%%
# =============================================================================
# redshift histograms
# =============================================================================
'''
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
'''

#%%
# =============================================================================
# bias test histograms
# =============================================================================
'''
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
'''



#%% 
# =============================================================================
# heatplots from GMM, 1.25<z<2 with different lower mass cuts
# =============================================================================
'''
num = 3
z_med_hp_low = 1.25
z_med_hp_high = 2.0
hp_lower_masses = [8.0, 8.5, 9.0]

z_med_hp = (z_med_hp_low+z_med_hp_high)/2.0
z_med_hp_gap = (z_med_hp_low+z_med_hp_high)/2.0 - z_med_hp_low
santini_idx = 1 # 1.3 to 2.0

for m in range(len(hp_lower_masses)):
    
    idx_rdm = np.arange(len(s))[(s['redshift_BEAGLE']>z_med_hp_low)&(s['redshift_BEAGLE']<z_med_hp_high)&(s['mass_BEAGLE_stellar']+s['mag_AD']>hp_lower_masses[m])] 
    print(len(idx_rdm))
    
    
    
    plt.hist((s['mass_BEAGLE_stellar']+s['mag_AD'])[idx_rdm], histtype='step', bins=30)
    plt.hist(s['mass_BEAGLE_stellar'][idx_rdm], histtype='step', bins=30)
    plt.show()
    
    
    
    plt.hist(s['sfr_BEAGLE_instant'][idx_rdm], histtype='step', bins=30)
    plt.show()    
    
    
    
    x_hp = np.array([])
    y_hp = np.array([])
    z_hp = np.array([])
    
    n_hp = 300 # number of samples to take from GMM in total
    
    for i in idx_rdm:

        for G in range(3):
            
            mean = np.array([s['x_GMM_3d'][i,G],s['y_GMM_3d'][i,G],s['z_GMM_3d'][i,G]])
            cov = np.array([[np.power(s['xsig_GMM_3d'][i,G],2), s['xycov_GMM_3d'][i,G], s['xzcov_GMM_3d'][i,G]],[s['xycov_GMM_3d'][i,G], np.power(s['ysig_GMM_3d'][i,G],2), s['yzcov_GMM_3d'][i,G]],[s['xzcov_GMM_3d'][i,G], s['yzcov_GMM_3d'][i,G], np.power(s['zsig_GMM_3d'][i,G],2)]])
    
            xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s['amp_GMM_3d'][i,G]))
    
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

    ax1.scatter(s['mass_BEAGLE_stellar'][idx_rdm], s['sfr_BEAGLE_instant'][idx_rdm], marker='x', color='r')


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
'''



#%%
'''
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
'''
#%%







