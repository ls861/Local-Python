#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 23:08:27 2021

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner

filenames = ['scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all', 
             'scenario_29_clusters_z2p0-3p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all', 
             'scenario_29_clusters_z3p0-4p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all',
             'scenario_29_clusters_z4p0-5p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all',
             'scenario_29_clusters_z5p0-6p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all',
             'scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm8p5',
             'scenario_29_clusters_z1p25-2p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm9p0']

#filenames = ['scenario_29_clusters_z1p25-2p0_4x32_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all']

filenames = ['scenario_29_clusters_z3p0-4p0_4x50000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_no_hogg']

filenames = ['scenario_29_clusters_z1p25-2p0_1x32_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_lm9real']

filenames = ['scenario_29_clusters_z1p25-2p0_4x5000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts', 
             'scenario_29_clusters_z1p25-2p0_4x5000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cutsNO']


filenames = ['scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts', 
             'scenario_29_clusters_z2p0-3p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts',
             'scenario_29_clusters_z3p0-4p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts',
             'scenario_29_clusters_z4p0-5p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts',
             'scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_lm8p5',
             'scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_lm9p0'
             ]

filenames = ['scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_limited_outlier', 
             'scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_limited_outlier_lm8p5', 
             'scenario_29_clusters_z1p25-2p0_4x15000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_limited_outlier_lm9p0']


filenames = ['lm_chain_3d_Hogg_truncated_normal_sSFR_scenario_29_clusters_z1p25-2p0_1x2000_true_linear_alpha_m_sfr_z_corr_alpha_beta_covprop_propscale_rdm_all_bad_gmm_cuts_limited_outlier_lm0p0.p']



filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_0_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_1_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_2_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_3_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_4_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_5_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_6_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_7_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_8_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x5000_9_001.p']

filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_12_001.p',
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x10000_12_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_12_001.p',
             'lm_chain_scenario_29_clusters_z4p0-5p0_4x10000_12_001.p',
             'lm_chain_scenario_29_clusters_z5p0-6p0_4x10000_12_001.p']


filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_14_001.p',
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x10000_14_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_14_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_14_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_16_001.p']

filenames = ['lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_12_001.p']

filenames = ['lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_17_001.p']




#%%
for filename in filenames:

    with open('/Users/LSand/Documents/linmix_files/{}'.format(filename), 'rb') as f:
     chain_original = pickle.load(f, encoding='latin1')
    

    # names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
#    names = ['pi', 'mu', 'tausqr']
    names = chain_original.dtype.names
    nChains=4
    minIter=len(chain_original)/nChains
    
#    burn=int(0.5*minIter)
    burn=0
    cap=minIter
    #cap=1900
    chain_arr = []
    for i in range(nChains):
        start = int(minIter*i+burn)
        finish = int((minIter*(i+1))-(minIter-cap))

        if filename == 'scenario_27_clusters_z1p25-2p0_4x5000':
            if i <= 1:
                chain_arr.append(chain_original[start:finish])
        else:
            if i < 20:
                chain_arr.append(chain_original[start:finish])                

    chain = np.concatenate(chain_arr)

    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        plt.plot(chain[name])
        plt.show()    
    print(' ')
    
#%%
'''
s = 0000
f = 3000

for i in range(len(chain['xi'][0])):


    plt.figure(figsize=(16, 4))
    plt.plot(chain['xi'][s:f,i])
    # plt.show()
    
    # plt.figure(figsize=(16, 4))
    plt.plot(chain['eta'][s:f,i])
    # plt.show()
    
    # plt.figure(figsize=(16, 4))
    plt.plot(chain['zeta'][s:f,i])
    plt.show()
    '''
#%%
for filename in filenames:
    with open('/Users/LSand/Documents/linmix_files/{}'.format(filename), 'rb') as f:
     chain_original = pickle.load(f, encoding='latin1')
    
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
#    names = ['pi', 'mu', 'tausqr']
#    names = chain_original.dtype.names
    nChains=4
    minIter=len(chain_original)/nChains
    
    burn=int(0.5*minIter)
#    burn=0
    cap=minIter
    #cap=1900
    chain_arr = []
    for i in range(nChains):
        start = int(minIter*i+burn)
        finish = int((minIter*(i+1))-(minIter-cap))

        if filename == 'scenario_27_clusters_z1p25-2p0_4x5000':
            if i <= 1:
                chain_arr.append(chain_original[start:finish])
        else:
            if i < 20:
                chain_arr.append(chain_original[start:finish])                

    chain = np.concatenate(chain_arr)

    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        plt.plot(chain[name])
        plt.show()    
    print(' ')
    pickle.dump(chain, open('/Users/LSand/Documents/linmix_files/PROCESSED_{}'.format(filename),'wb'))


    idx = np.random.choice(len(chain['alphaN_a']), size=1000, replace=True)
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']

    sample_chain = {}
    for name in names:
        sample_chain[name] = chain[name][idx]
    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
        plt.plot(sample_chain[name])
        plt.show()
    pickle.dump(sample_chain, open('/Users/LSand/Documents/linmix_files/1000_sample_chains/PROCESSED_{}'.format(filename),'wb'))


    #%%
    ## =============================================================================
    ## Corner Plots
    ## =============================================================================

    chain_corner = sample_chain
    
    names_plot = ['ssfr a', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
    data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['beta_a'][:], chain_corner['sig0'][:], chain_corner['pbad'][:], chain_corner['outlier_mean'][:], chain_corner['outlier_sigma'][:]]).T

#    names_plot = ['ssfr a', 'beta a', 'sig0', 'outlier mean', 'outlier sigma']
#    data_corner = np.array([chain_corner['alphaN_a'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
#
#    names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
#    data_corner = np.array([chain_corner['alphaN_a'], chain_corner['alphaN_b'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['pbad'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
#
#    names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'outlier mean', 'outlier sigma']
#    data_corner = np.array([chain_corner['alphaN_a'], chain_corner['alphaN_b'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
#    
#    names_plot = ['ssfr a', 'beta a', 'sig0', 'k', 'outlier mean', 'outlier sigma']
#    data_corner = np.array([chain_corner['alphaN_a'], chain_corner['beta_a'], chain_corner['sig0'], chain_corner['k'], chain_corner['outlier_mean'], chain_corner['outlier_sigma']]).T
        
    figure = corner.corner(data_corner, labels=names_plot,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True, title_kwargs={"fontsize": 20})



