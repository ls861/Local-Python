#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 11:07:14 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# =============================================================================
# choices
# =============================================================================

# never made mass or sfr inputs for mock, hence NA for MS options

option_mock                         = 1 # 1 is usual mock, 2 inc GMMz, GMMxsig, GMMxzcov, GMMyzcov
option_main_sequence                = 0 # NA for mock
option_chains                       = 1 # 1, 2 is redshift dependence hack
option_chains_2x2                   = 0
option_chains_4x1                   = 0
option_ind_objects_mass             = 1
option_ind_objects_sfr              = 1
option_ind_objects_z                = 0 # 2 only, redshift dependence hack
option_samples_from_fitted_model    = 0
option_fitted_model_over_ms         = 0 # NA for mock
option_heatplots                    = 0 # legend labels not auto update


# =============================================================================
# load inputs
# =============================================================================

# if I want to view the mock results - JUST chains
# 002 4 500
# 001 4 5000
# 001 4 50000

#if I want ASTRODEEP, first 4 fields
# 8p5 2p0 0p5 9p5 002 removed mass > 10.5, SFR < -2 and SFR > 2.1, first 4 fields      
# 002 4 500
# 001 4 5000

# if I want to view the mock results - JUST chains for redshift dependent alpha beta runs
# 001 1 20 
# 001 1 200
# 001 1 2000
# 001 4 20
# 001 4 200 
# 001 4 2000 

# actually added in alpha beta redshift dependence...
# 001 4 30
# 001 4 300 
# 001 4 3000 

# 002 4 3000 (2) NO mass dependence of scatter, or z dependence of alpha and beta

# 003 4 600 (1) checking 2d script removal of redshift dependence works

comment2    = '003'
nChains     = 4             # options are 1 or 4
minIter     = 600            # options are 20, 200, 2000, or 
burn        = 0



if option_mock == 1:
    sbf = '/Users/lester/Documents/GitHub/M_SFR_project_routines/limix_inputs_106_x47_mass5p0/' # subfolder for DE 106, 47 objects, mass > 5 (ie all of them)
    
if option_mock == 2:
    sbf = '/Users/lester/Documents/GitHub/M_SFR_project_routines/linmix_inputs_z_dependence/7p0_' # subfolder for DE 106, 47 objects, mass > 5 (ie all of them), includes GMMz, GMMxsig, GMMxzcov, GMMyzcov
    
if option_mock == 1 or option_mock == 2:

    pi_err      = np.load(sbf+'pi_err.npy')         # 3x probability of each posterior gaussian
    GMMx        = np.load(sbf+'GMMx.npy')           # 3x posterior means per mass
    GMMy        = np.load(sbf+'GMMy.npy')           # 3x posterior means per sfr
    GMMxsig     = np.load(sbf+'GMMxsig.npy')        # 3x posterior sigmas per mass
    GMMysig     = np.load(sbf+'GMMysig.npy')        # 3x posterior sigmas per sfr
    GMMxycov    = np.load(sbf+'GMMxycov.npy')       # 3x posterior covar per mass-sfr pair
    
    nK          = np.load(sbf+'nK.npy')             # 3 #gaussians modelling xi
    nGauss      = np.load(sbf+'nGauss.npy')         # 3 #gaussians modelling BEAGLE posterior
    
#    nChains     = np.load(sbf+'nChains.npy')        # 2
#    minIter     = np.load(sbf+'minIter.npy')        # 3000
#    maxIter     = np.load(sbf+'maxIter.npy')        # 3000
 
# NOT CREATED FOR MOCK YET
#    GMMz        = np.load(sbf+'GMMz.npy')  
#    GMMchi2     = np.load(sbf+'GMMchi2.npy')  
#    GMMmass     = np.load(sbf+'GMMmass.npy')  
#    GMMsfr      = np.load(sbf+'GMMsfr.npy')     
    
    # TRUE VALUES OF RELATION
    true_values = {'alpha':-6.2, 'beta':0.8, 'sig_0':0.3, 'k':1.0}
    
elif option_mock == 0:
    #massLim = '8p5'
    #zLim = '2p0'
    #wLim = '0p5'
    #chi2Lim = '9p5'
    #
    #massLim = '7p0'
    #zLim = '7p0'
    #wLim = '10p0'
    #chi2Lim = '9999'
    
    massLim = '8p5'
    zLim = '2p0'
    wLim = '0p5'
    chi2Lim = '9p5'
    comment = '002'
    
    # 8p5 2p0 0p5 9p5 002 removed mass > 10.5, SFR < -2 and SFR > 2.1, first 4 fields      
    sbf = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/linmix_npy_files_z{}_w{}_chi{}_{}/{}_'.format(zLim, wLim, chi2Lim, comment, massLim)
    
    pi_err      = np.load(sbf+'pi_err.npy')         # 3x probability of each posterior gaussian
    GMMx        = np.load(sbf+'GMMx.npy')           # 3x posterior means per mass
    GMMy        = np.load(sbf+'GMMy.npy')           # 3x posterior means per sfr
    GMMxsig     = np.load(sbf+'GMMxsig.npy')        # 3x posterior sigmas per mass
    GMMysig     = np.load(sbf+'GMMysig.npy')        # 3x posterior sigmas per sfr
    GMMxycov    = np.load(sbf+'GMMxycov.npy')       # 3x posterior covar per mass-sfr pair
    
    nK          = np.load(sbf+'nK.npy')             # 3 #gaussians modelling xi
    nGauss      = np.load(sbf+'nGauss.npy')         # 3 #gaussians modelling BEAGLE posterior
    
#    nChains     = np.load(sbf+'nChains.npy')        # 2
#    minIter     = np.load(sbf+'minIter.npy')        # 3000
#    maxIter     = np.load(sbf+'maxIter.npy')        # 3000
    
    GMMz        = np.load(sbf+'GMMz.npy')  
    GMMchi2     = np.load(sbf+'GMMchi2.npy')  
    GMMmass     = np.load(sbf+'GMMmass.npy')  
    GMMsfr      = np.load(sbf+'GMMsfr.npy')  




# =============================================================================
# MAIN SEQUENCE
# =============================================================================

if option_main_sequence == 1 and option_mock != 1:
 
    print('MAIN SEQUENCE')
    
    plt.xlim(7, 11)
    plt.ylim(-2, 3)
    plt.scatter(GMMmass, GMMsfr, marker='x')
    plt.show()
    
    plt.xlim(7, 12)
    plt.ylim(-5, 5)
    plt.scatter(GMMmass, GMMsfr, marker='x')
    plt.show()
    
    plt.scatter(GMMmass, GMMsfr, marker='x')
    plt.xlabel(r'$(m_{tot}/M_{\odot})$')
    
    plt.show()

# =============================================================================
# load finished chain
# =============================================================================

#lm = np.load('lm_30000.npy')[5000:]
#lm = np.load('/Users/lester/Documents/GitHub/M_SFR_project_routines/linmix_inputs_TEST/7p5_lm_chain.npy')
#lm = np.load('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/linmix_npy_files_z2p0_chi9p5/8p5_lm_chain_30000.npy')
#lm = np.load('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jun_2020/linmix_npy_files_z2p0_chi9p5/8p5_lm_chain_4x_10000_001.npy')
#lm = np.load(sbf+'lm_chain_4x_10000_001.npy')
    
if option_mock == 2:
    lm = np.load(sbf+'lm_chain_{}x_{}_{}_z_dependence.npy'.format(str(nChains), str(minIter), comment2))
else:
    lm = np.load(sbf+'lm_chain_{}x_{}_{}.npy'.format(str(nChains), str(minIter), comment2))
    

# =============================================================================
# combine chains and remove burn in
# =============================================================================

lm_arr = []
for i in range(nChains):
    start = minIter*i+burn
    finish = minIter*(i+1)
    lm_arr.append(lm[start:finish])
lm = np.concatenate(lm_arr)

# =============================================================================
# PLOTS of alpha beta sig_0 k
# =============================================================================

params = ['alpha', 'beta', 'sig_0', 'k']

if option_chains == 1:

    print('ALPHA BETA SIG0 K')
    
    for i, param in enumerate(params):
    
        plt.title(param.replace('_', ''))
        plt.plot(lm[param])
        plt.show()
        
        plt.figure(figsize=(6, 6))
    #    plt.title('Histogram representing the posterior distribution of {}'.format(param.replace('_', '')))
        y, x, _ = plt.hist(lm[param], bins=20)
        plt.xlabel(param.replace('_', ''))
        plt.ylabel('Count')
        plt.plot((np.median(lm[param]), np.median(lm[param])),(0, y.max()), color='r', label='median')
        if option_mock == 1:
            plt.plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        plt.legend()
        plt.show()


# =============================================================================
# PLOTS of alpha beta sig_0 k - REDSHIFT DEPENDENCE TEST
# =============================================================================

if option_chains == 2:
    print('ALPHA BETA SIG0 K')
    # takes first object only for alpha and beta chain
    
    idx_i = 0
    
    params = ['alpha', 'beta']
    for i, param in enumerate(params):
    
        plt.title(param.replace('_', ''))
        plt.plot(lm[param][:,idx_i])
        plt.show()
        
        plt.figure(figsize=(6, 6))
    #    plt.title('Histogram representing the posterior distribution of {}'.format(param.replace('_', '')))
        y, x, _ = plt.hist(lm[param][:,idx_i], bins=20)
        plt.xlabel(param.replace('_', ''))
        plt.ylabel('Count')
        plt.plot((np.median(lm[param][:,idx_i]), np.median(lm[param][:,idx_i])),(0, y.max()), color='r', label='median')
        if option_mock == 1 or option_mock == 2:
            plt.plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        plt.legend()
        plt.show()
    
    params = ['sig_0', 'k']
    for i, param in enumerate(params):
    
        plt.title(param.replace('_', ''))
        plt.plot(lm[param])
        plt.show()
        
        plt.figure(figsize=(6, 6))
    #    plt.title('Histogram representing the posterior distribution of {}'.format(param.replace('_', '')))
        y, x, _ = plt.hist(lm[param], bins=20)
        plt.xlabel(param.replace('_', ''))
        plt.ylabel('Count')
        plt.plot((np.median(lm[param]), np.median(lm[param])),(0, y.max()), color='r', label='median')
        if option_mock == 1 or option_mock == 2:
            plt.plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        plt.legend()
        plt.show()
        
    params = ['alpha_a', 'alpha_b', 'alpha_c', 'beta_a', 'beta_b', 'beta_c']
    plt.figure(figsize=(10, 10))
    for i, param in enumerate(params):
    
#        plt.title(param.replace('_',''))
        plt.plot(lm[param], label=param.replace('_',''))
    
    plt.legend()
    plt.show()
        
    for i, param in enumerate(params):
        
        plt.figure(figsize=(6, 6))
    #    plt.title('Histogram representing the posterior distribution of {}'.format(param.replace('_', '')))
        y, x, _ = plt.hist(lm[param], bins=20)
        plt.xlabel(param.replace('_', ''))
        plt.ylabel('Count')
        plt.plot((np.median(lm[param]), np.median(lm[param])),(0, y.max()), color='r', label='median')
#        if option_mock == 1 or option_mock == 2:
#            plt.plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        plt.legend()
        plt.show()


# =============================================================================
# histograms but in 2x2
# =============================================================================

if option_chains_2x2 == 1:

    fig, axs = plt.subplots(2, 2, figsize=(14,10))
    #fig.suptitle('Posterior distributions for each fitted parameter', size=20)
    #plt.figure(figsize=(10, 10))
    
    k = 0
    j = 0    
    
    for i, param in enumerate(params):
    
        y, x, _ = axs[k,j].hist(lm[param], bins=20)
        axs[k,j].set_xlabel(param.replace('_', ''))
        axs[k,j].set_ylabel('Count')
        axs[k,j].plot((np.median(lm[param]), np.median(lm[param])),(0, y.max()), color='r', label='median')
        if option_mock == 1:
            axs[k,j].plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        axs[k,j].legend()
    
        
        if k==0 and j==0:
            j=1
        elif k==0 and j==1:
            k=1
            j=0
        elif k==1 and j==0:
            j=1
    
    plt.show()

# =============================================================================
# histograms but in 4x1
# =============================================================================


if option_chains_4x1 == 1:

    fig, axs = plt.subplots(1, 4, figsize=(24,3))
    #fig.suptitle('Posterior distributions for each fitted parameter', size=20)
    #plt.figure(figsize=(10, 10))
    
    xlabels = [r'$\alpha$', r'$\beta$', r'$\sigma_0$', r'$\kappa$',]
    
    k = 0
    j = 0   
    
    for i, param in enumerate(params):
    
        y, x, _ = axs[k].hist(lm[param], bins=20)
    #    axs[k].set_xlabel(param.replace('_', ''))
        axs[k].set_xlabel(xlabels[i])
        axs[k].set_ylabel('Count')
        axs[k].plot((np.median(lm[param]), np.median(lm[param])),(0, y.max()), color='r', label='median')
        if option_mock == 1:
            axs[k].plot((true_values[param], true_values[param]),(0, y.max()), color='k', label='true')
        axs[k].legend()
    
        k+=1
    
    plt.show()


# =============================================================================
# plotting xi chain for sample of objects
# =============================================================================

if option_ind_objects_mass == 1:
    
    plt.figure(figsize=(8, 5))
    for i in range(len(lm['xi'][0,21:33])):
        plt.plot(lm['xi'][:,i])
    plt.show()

# =============================================================================
# plotting eta chain for sample of objects
# =============================================================================

if option_ind_objects_sfr == 1:
    
    plt.figure(figsize=(8, 5))
    for i in range(len(lm['eta'][0,21:33])):
        plt.plot(lm['eta'][:,i])
    plt.show()
    
# =============================================================================
# plotting zeta chain for sample of objects
# =============================================================================

if option_ind_objects_z == 2:
    
    plt.figure(figsize=(8, 5))
    for i in range(len(lm['zeta'][0,21:33])):
        plt.plot(lm['zeta'][:,i])
    plt.show()
    
# =============================================================================
# Fitted model over Main Sequence
# =============================================================================

xlow = 8.5
xhigh = 10.0
count = 10000
x = np.linspace(xlow, xhigh, count)

alpha = np.median(lm['alpha'])
beta = np.median(lm['beta'])
sig_0 = np.median(lm['sig_0'])
k = np.median(lm['k'])
    
if option_samples_from_fitted_model == 1:

    print('SAMPLES TAKEN FROM FITTED MODEL')
    print(alpha, beta, sig_0, k)
    
    
    xi_min = 8.5
    xi_max = 10.0
    
    sigma = sig_0 * ( ((1.0-k)*(x-xi_max) / (xi_max - xi_min)) + 1)
    
    y = alpha + (beta*x) + np.random.normal(0, abs(sigma), len(x))
    
    plt.hist2d(x, y, bins=50)
    plt.show()
    
    plt.xlim(xlow, xhigh)
    plt.ylim(-1, 4)
    plt.scatter(x, y, marker='x')
    plt.show()

# =============================================================================
# plotting fitted params on input main sequence
# =============================================================================

if option_fitted_model_over_ms and option_mock != 1:

    print('FITTED MODEL OVER MAIN SEQUENCE')
    
    # SANTINI, log(SFR) = alpha log(M / M_9p7) + beta
    # 1.3 < z < 2
    s1_alpha = 1.04
    s1_beta = 1.01
    s1_alpha_err = 0.03
    s1_beta_err = 0.04
    # 2.0 < z < 3.0
    s2_alpha = 1.16
    s2_beta = 1.22
    s2_alpha_err = 0.03
    s2_beta_err = 0.03
    
    xlow = 8.0
    xhigh = 10.5
    ylow = -1.5
    yhigh = 2
    siglow = sig_0 * ( ((1.0-k)*(xlow-xi_max) / (xi_max - xi_min)) + 1)
    sighigh = sig_0 * ( ((1.0-k)*(xhigh-xi_max) / (xi_max - xi_min)) + 1)
    
    plt.figure(figsize=(10, 6))
    plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
    plt.xlim(xlow, xhigh)
    plt.ylim(ylow, yhigh)
    plt.scatter(GMMmass, GMMsfr, marker='x')
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), color='k')
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='r')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='r')
    
    # santini z=1
    s_x = np.array([8.4, 9.2, 10.0])
    s_y = s1_alpha*(s_x - 9.7) + s1_beta
    s_y_intrinsic = np.array([0.36, 0.35, 0.0])
    s_y_observed = np.array([0.51, 0.46, 0.26])
    
    print('alpha', s1_alpha*(0 - 9.7) + s1_beta)
    print('beta', (s1_alpha*(1 - 9.7) + s1_beta) - (s1_alpha*(0 - 9.7) + s1_beta) )
    
    plt.plot(s_x, s_y, color='g', label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    plt.plot(s_x, s_y+s_y_intrinsic, color='g', linestyle='--', label=r'Santini+17, Intrinsic Scatter')
    plt.plot(s_x, s_y-s_y_intrinsic, color='g', linestyle='--')
    plt.plot(s_x, s_y+s_y_observed, color='g', linestyle=':', label=r'Santini+17, Observed Scatter')
    plt.plot(s_x, s_y-s_y_observed, color='g', linestyle=':')
    
    # santini z=1 WRONG
    #plt.plot((xlow, xhigh), (s1_alpha*(xlow - 9.7) + s1_beta, s1_alpha*(xhigh - 9.7) + s1_beta), color='g', label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    #s_low1 = (s1_alpha+s1_alpha_err)*(xlow - 9.7) + (s1_beta+s1_beta_err)
    #s_low2 = (s1_alpha-s1_alpha_err)*(xlow - 9.7) + (s1_beta-s1_beta_err)
    #s_low3 = (s1_alpha+s1_alpha_err)*(xlow - 9.7) + (s1_beta-s1_beta_err)
    #s_low4 = (s1_alpha-s1_alpha_err)*(xlow - 9.7) + (s1_beta+s1_beta_err)
    #s_high1 = (s1_alpha+s1_alpha_err)*(xhigh - 9.7) + (s1_beta+s1_beta_err)
    #s_high2 = (s1_alpha-s1_alpha_err)*(xhigh - 9.7) + (s1_beta-s1_beta_err)
    #s_high3 = (s1_alpha+s1_alpha_err)*(xhigh - 9.7) + (s1_beta-s1_beta_err)
    #s_high4 = (s1_alpha-s1_alpha_err)*(xhigh - 9.7) + (s1_beta+s1_beta_err)
    #s_low_min = min(s_low1, s_low2, s_low3, s_low4)
    #s_low_max = max(s_low1, s_low2, s_low3, s_low4)
    #s_high_min = min(s_high1, s_high2, s_high3, s_high4)
    #s_high_max = max(s_high1, s_high2, s_high3, s_high4)
    #plt.plot((xlow, xhigh), (s_low_max, s_high_max), color='g', linestyle='--')
    #plt.plot((xlow, xhigh), (s_low_min, s_high_min), color='g', linestyle='-.')
    
    # santini z=2 WRONG
    #plt.plot((xlow, xhigh), (s2_alpha*(xlow - 9.7) + s2_beta, s2_alpha*(xhigh - 9.7) + s2_beta), color='b', label=r'Santini+17, 2.0 $<$ z $<$ 3.0')
    #s_low1 = (s2_alpha+s2_alpha_err)*(xlow - 9.7) + (s2_beta+s2_beta_err)
    #s_low2 = (s2_alpha-s2_alpha_err)*(xlow - 9.7) + (s2_beta-s2_beta_err)
    #s_low3 = (s2_alpha+s2_alpha_err)*(xlow - 9.7) + (s2_beta-s2_beta_err)
    #s_low4 = (s2_alpha-s2_alpha_err)*(xlow - 9.7) + (s2_beta+s2_beta_err)
    #s_high1 = (s2_alpha+s2_alpha_err)*(xhigh - 9.7) + (s2_beta+s2_beta_err)
    #s_high2 = (s2_alpha-s2_alpha_err)*(xhigh - 9.7) + (s2_beta-s2_beta_err)
    #s_high3 = (s2_alpha+s2_alpha_err)*(xhigh - 9.7) + (s2_beta-s2_beta_err)
    #s_high4 = (s2_alpha-s2_alpha_err)*(xhigh - 9.7) + (s2_beta+s2_beta_err)
    #s_low_min = min(s_low1, s_low2, s_low3, s_low4)
    #s_low_max = max(s_low1, s_low2, s_low3, s_low4)
    #s_high_min = min(s_high1, s_high2, s_high3, s_high4)
    #s_high_max = max(s_high1, s_high2, s_high3, s_high4)
    #plt.plot((xlow, xhigh), (s_low_max, s_high_max), color='b', linestyle='--')
    #plt.plot((xlow, xhigh), (s_low_min, s_high_min), color='b', linestyle='--')
    
    plt.legend()
    plt.show()

# =============================================================================
# heatplots
# =============================================================================

if option_heatplots == 1:
    
    print('HEATPLOTS')
    
    hp_x = []
    hp_y = []
    
    samples = 100
    
    # objs to plot
    #objs = np.random.choice(range(len(GMMx)), 2) # random selection
    objs = range(len(GMMx)) # all
    
    for obj in objs:
    #    print(obj)
        draws = np.random.choice([0, 1, 2], samples, p=pi_err[obj]/sum(pi_err[obj]))
        for draw in draws:
            mean = (GMMx[obj,draw], GMMy[obj,draw])
            cov = ((GMMxsig[obj,draw],GMMxycov[obj,draw]),(GMMxycov[obj,draw],GMMysig[obj,draw]))
            hpx, hpy = np.random.multivariate_normal(mean, cov)
            hp_x.append(hpx)
            hp_y.append(hpy)
    
    
    plt.figure(figsize=(10, 6))
    plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
    plt.hist2d(hp_x, hp_y, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
    plt.colorbar()
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax = plt.axes()
    ax.set_facecolor(rgba)
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), color='k')
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='r')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='r')
    plt.show()
    
    
    # DECENT PLOT
    
    plt.figure(figsize=(10, 6))
    lw = 4
    #plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
    plt.hist2d(hp_x, hp_y, bins=50, range=((xlow, xhigh),(ylow, yhigh)), cmap=plt.cm.viridis)
    plt.colorbar(label='Count')
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax = plt.axes()
    ax.set_facecolor(rgba)
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), linewidth=lw, color='k', label=r'Fitted, 1.5 $<$ z $<$ 2.5')
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), linestyle='--', linewidth=lw, color='k', label=r'Fitted, Intrinsic Scatter')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), linestyle='--', linewidth=lw, color='k')
    
    plt.plot(s_x, s_y, color='w', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    plt.plot(s_x, s_y+s_y_intrinsic, color='w', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    plt.plot(s_x, s_y-s_y_intrinsic, color='w', linestyle='--', linewidth=lw)
    plt.plot(s_x, s_y+s_y_observed, color='w', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    plt.plot(s_x, s_y-s_y_observed, color='w', linestyle=':', linewidth=lw)
    
    plt.legend(loc='lower right')
    plt.show()    
    
    plt.scatter(hp_x, hp_y)
    plt.show()
    
    












