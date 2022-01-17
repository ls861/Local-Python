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


filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x20000_18_001.p']



filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_18_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_19_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_20_001.p',
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_21_001.p',

             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_18_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_19_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_20_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_21_001.p',
             
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_20_001.p',
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_21_001.p',
             'lm_chain_scenario_29_clusters_z4p0-5p0_4x50000_21_001.p']


filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_21_001.p',
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_21_001.p',
             'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_21_001.p',
             'lm_chain_scenario_29_clusters_z4p0-5p0_4x50000_21_001.p']


filenames = ['lm_chain_scenario_31_clusters_z1p25-2p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_21_001.p']

filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_19_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_19_001.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_20_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p']


filenames = ['lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_19_001.p']




filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_19_001.p',
             'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_19_001.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_19_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_19_001.p',
             
             'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_20_001.p',
             'lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_20_001.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_20_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p',

             'lm_chain_scenario_31_clusters_z1p25-2p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_001.p',
             'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_21_001.p']


# filenames = ['lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p',

#              'lm_chain_scenario_31_clusters_z1p25-2p0_4x50000_21_001.p',
#              'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_21_001.p',
#              'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_001.p',
#              'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_21_001.p']


filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_002.p',
             'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_003.p']

filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_004.p']


filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_23_001.p']
# filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_20_001.p']
# filenames = ['lm_chain_scenario_31_clusters_z3p0-4p0_1x5000_23_001.p']

filenames = ['lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_22_001.p']

filenames = ['lm_chain_scenario_34_clusters_z1p25-6p0_4x50_24_001.p']
filenames = ['lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_24_001.p']

filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x50_27_001.p']

filenames = ['lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_22_001.p',
             'lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_24_001.p',
             'lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_24_002.p',
             'lm_chain_scenario_34_clusters_z1p25-6p0_4x5000_25_001.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x5000_27_001.p']

filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x5000_26_001.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x5000_26_002.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x5000_26_003.p']


filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x10000_27_001.p']


filenames = ['lm_chain_scenario_34_clusters_z1p25-5p0_4x5000_24_003.p',
             'lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_004.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x10000_28_001.p']


filenames = ['lm_chain_scenario_34_clusters_z3p0-5p0_4x5000_24_003.p']

filenames = ['lm_chain_scenario_34_clusters_z1p25-3p0_4x5000_24_003.p',
             'lm_chain_scenario_34_clusters_z2p0-4p0_4x5000_24_003.p',
             'lm_chain_scenario_34_clusters_z3p0-5p0_4x5000_24_003.p',
             'lm_chain_scenario_34_clusters_z1p25-4p0_4x5000_24_003.p',
             'lm_chain_scenario_34_clusters_z2p0-5p0_4x5000_24_003.p']

filenames = ['lm_chain_scenario_34_clusters_z1p25-3p0_4x5000_24_003.p',
             'lm_chain_scenario_34_clusters_z3p0-5p0_4x5000_24_003.p']


filenames = ['lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_004.p']

filenames = ['lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_20_001.p']
filenames = ['lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p']


# filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_26_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_26_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_26_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_26_005.p']
filenames = ['lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_005.p']

# filenames = ['lm_chain_scenario_33_clusters_z1p25-6p0_4x20000_24_005.p']
# filenames = ['lm_chain_scenario_34_clusters_z1p25-6p0_4x20000_24_005.p']


# PBAD = 0
# filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_28_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_28_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_28_005.p']
# filenames = ['lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_28_005.p']
filenames = ['lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_28_005.p']


# PBAD = 0, + sigma clipping
# filenames = ['lm_chain_scenario_35_clusters_z1p25-2p0_4x20000_28_006.p']
# filenames = ['lm_chain_scenario_35_clusters_z2p0-3p0_4x20000_28_006.p']
# filenames = ['lm_chain_scenario_35_clusters_z3p0-4p0_4x20000_28_006.p']
# filenames = ['lm_chain_scenario_35_clusters_z4p0-5p0_4x20000_28_006.p']
# filenames = ['lm_chain_scenario_35_clusters_z5p0-6p0_4x20000_28_006.p']


filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_28_005.p',
             'lm_chain_scenario_35_clusters_z1p25-2p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z2p0-3p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z3p0-4p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z4p0-5p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z5p0-6p0_4x20000_28_006.p',
             'lm_chain_scenario_33_clusters_z1p25-6p0_4x20000_24_005.p',
             'lm_chain_scenario_34_clusters_z1p25-6p0_4x20000_24_005.p']
             
             



# FULL CHAINS

#%%
for filename in filenames:

    with open('/Users/lester/Documents/linmix_files/{}'.format(filename), 'rb') as f:
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
    with open('/Users/lester/Documents/linmix_files/{}'.format(filename), 'rb') as f:
        chain_original = pickle.load(f, encoding='latin1')
    
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']
    # names = ['alphaN_a', 'beta_a', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
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
                
        elif filename == 'lm_chain_scenario_34_clusters_z1p25-3p0_4x5000_24_003.p':
            if i <= 2:
                chain_arr.append(chain_original[start:finish])
                
        elif filename == 'lm_chain_scenario_34_clusters_z3p0-5p0_4x5000_24_003.p':
            if i == 3:
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
    
    # pickle.dump(chain, open('/Users/lester/Documents/linmix_files/PROCESSED_{}'.format(filename),'wb'))

    # idx = np.random.choice(len(chain['alphaN_a']), size=len(chain['alphaN_a']), replace=False)
    names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'k', 'pbad', 'outlier_mean', 'outlier_sigma']

    sample_chain = {}
    for name in names:
        sample_chain[name] = chain[name]
    # for name in names:            
    #     plt.figure(figsize=(16, 4))
    #     plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
    #     plt.plot(sample_chain[name])
    #     plt.show()
    # pickle.dump(sample_chain, open('/Users/lester/Documents/linmix_files/simple_chains/PROCESSED_{}'.format(filename),'wb'))


    #%%
    
    
    ## =============================================================================
    ## Corner Plots
    ## =============================================================================

for filename in filenames:
    # with open('/Users/lester/Documents/linmix_files/100000_sample_chains/PROCESSED_{}'.format(filename), 'rb') as f:
    #     chain_corner = pickle.load(f, encoding='latin1')
    with open('/Users/lester/Documents/linmix_files/PROCESSED_{}'.format(filename), 'rb') as f:
        chain_corner = pickle.load(f, encoding='latin1')    

    if '_19_' in filename or '_21_' in filename or '_28_' in filename:
    
        # filename == 'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_19_001.p' or \
        # filename == 'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_19_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_19_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_19_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_19_001.p' or \
        # filename == 'lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_29_clusters_z4p0-5p0_4x50000_21_001.p'or \
        # filename == 'lm_chain_scenario_31_clusters_z1p25-2p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z2p0-3p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_21_001.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_002.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_003.p' or \
        # filename == 'lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_21_004.p':     
            
        names_plot = ['ssfr a', 'beta a', 'sig0']
        data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['beta_a'][:], chain_corner['sig0'][:]]).T
    
    elif '_22_' in filename:
        names_plot = ['ssfr a', 'ssfr b', 'beta a', 'beta_b', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['alphaN_b'][:], chain_corner['beta_a'][:], chain_corner['beta_b'][:], chain_corner['sig0'][:], chain_corner['pbad'][:], chain_corner['outlier_mean'][:], chain_corner['outlier_sigma'][:]]).T

    elif '_24_' in filename or '_25_' in filename:
        names_plot = ['ssfr a', 'ssfr b', 'beta a', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['alphaN_b'][:], chain_corner['beta_a'][:], chain_corner['sig0'][:], chain_corner['pbad'][:], chain_corner['outlier_mean'][:], chain_corner['outlier_sigma'][:]]).T
        
    elif '_27_' in filename:
        names_plot = ['ssfr a', 'beta a', 'sig0', 'k', 'pbad', 'outlier mean', 'outlier sigma']
        data_corner = np.array([chain_corner['alphaN_a'][:], chain_corner['beta_a'][:], chain_corner['sig0'][:], chain_corner['k'][:], chain_corner['pbad'][:], chain_corner['outlier_mean'][:], chain_corner['outlier_sigma'][:]]).T
        
    else:
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


    plt.suptitle('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '), fontsize=20, y=1.1)
    plt.show()
    
    
    
    
    
    
    
    
#%%
    
# =============================================================================
# converged mass sfr and redshift chains
# =============================================================================
    
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner


filenames = ['lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_26_005.p',
             'lm_chain_scenario_33_clusters_z1p25-2p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z2p0-3p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z3p0-4p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_28_005.p',
             'lm_chain_scenario_33_clusters_z5p0-6p0_4x20000_28_005.p',
             'lm_chain_scenario_35_clusters_z1p25-2p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z2p0-3p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z3p0-4p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z4p0-5p0_4x20000_28_006.p',
             'lm_chain_scenario_35_clusters_z5p0-6p0_4x20000_28_006.p',
             'lm_chain_scenario_33_clusters_z1p25-6p0_4x20000_24_005.p',
             'lm_chain_scenario_34_clusters_z1p25-6p0_4x20000_24_005.p']


filenames = ['lm_chain_scenario_33_clusters_z4p0-5p0_4x20000_26_005.p']

for filename in filenames:
    with open("/Users/lester/Documents/linmix_files/PROCESSED_{}".format(filename), 'rb') as f:
        processed = pickle.load(f, encoding='latin1')          
    
    print(processed.dtype.names)
    
    
    names = ['xi', 'eta', 'zeta']
    
    sample_chain = {}
    for name in names:
        sample_chain[name] = processed[name]
        
    for name in names:            
        plt.figure(figsize=(16, 4))
        plt.title('{}\n{}'.format(filename, name).replace('_',' '))
        plt.plot(sample_chain[name])
        plt.show()


    arr_means = {}
    arr_medians = {}
    arr_medians5000 = {}

    for name in names:
        means = []
        medians = []
        medians5000 = []
        for i in range(len(sample_chain['xi'][0])):   
            # plt.figure(figsize=(16, 4))
            # plt.title('{}\n{} {} {:.2f} {:.2f}'.format(filename, name, i, np.mean(sample_chain[name][:,i]), np.median(sample_chain[name][:,i])).replace('_',' '))
            # plt.plot(sample_chain[name][:,i])
            # plt.show()       
            
            means.append(np.mean(sample_chain[name][:,i]))
            medians.append(np.median(sample_chain[name][:,i]))
            medians5000.append(np.median(sample_chain[name][:5000,i]))
        
        plt.scatter(means, medians)
        plt.scatter(means, medians5000, marker='x')
        plt.title('{}\n{}'.format(filename, name))
        plt.show()
        
        
        arr_means[name] = means
        arr_medians[name] = medians
        arr_medians5000[name] = medians5000   
        
        
    pickle.dump(arr_means, open('/Users/lester/Documents/linmix_files/mass_sfr_z_chains/means_{}'.format(filename),'wb'))
    pickle.dump(arr_medians, open('/Users/lester/Documents/linmix_files/mass_sfr_z_chains/medians_{}'.format(filename),'wb'))
    pickle.dump(arr_medians5000, open('/Users/lester/Documents/linmix_files/mass_sfr_z_chains/medians5000_{}'.format(filename),'wb'))


#%%

filename = filenames[0]
with open('/Users/lester/Documents/linmix_files/mass_sfr_z_chains/medians_{}'.format(filename), 'rb') as f:
    tt = pickle.load(f, encoding='latin1')   


print(tt['xi'])

plt.scatter(tt['xi'], tt['eta'])





#%%

# =============================================================================
# 
# =============================================================================

# FULL CHAINS SINGLE PLOTS
'''
#%%
for filename in filenames:

    with open('/Users/lester/Documents/linmix_files/{}'.format(filename), 'rb') as f:
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
        if len(np.shape(chain[name])) == 2:
            for n in range(np.shape(chain[name])[1]):
                
                plt.figure(figsize=(16, 4))
                plt.title('{}\n{} {}, burn:{}, cap:{}'.format(filename, name, n, burn, cap).replace('_',' '))
                plt.plot(chain[name][:,n])
                plt.show()    
                print(' ')
            
        else:
            plt.figure(figsize=(16, 4))
            plt.title('{}\n{}, burn:{}, cap:{}'.format(filename, name, burn, cap).replace('_',' '))
            plt.plot(chain[name])
            plt.show()    
            print(' ')
    
# print(names)
'''

    
'''
idx_temp = chain['tausqr'][:,2] > 10
print(chain['pi'][:,1][idx_temp])
print(chain['mu'][:,1][idx_temp])
print(chain['tausqr'][:,1][idx_temp])


plt.hist(chain['pi'][:,1][idx_temp], bins = 50)




filenames = ['PROCESSED_lm_chain_scenario_29_clusters_z1p25-2p0_4x50000_20_001.p',
             'PROCESSED_lm_chain_scenario_29_clusters_z2p0-3p0_4x50000_20_001.p',
             'PROCESSED_lm_chain_scenario_31_clusters_z3p0-4p0_4x50000_20_001.p',
             'PROCESSED_lm_chain_scenario_31_clusters_z4p0-5p0_4x50000_20_001.p']


for filename in filenames:

    with open('/Users/lester/Documents/linmix_files/{}'.format(filename), 'rb') as f:
        chain_original = pickle.load(f, encoding='latin1')
    
    print(np.shape(chain_original['xi']))

'''







    