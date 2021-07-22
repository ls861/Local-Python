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

option_main_sequence                = 1 # NA for mock
option_chains                       = 1 # 1, 2 is redshift dependence hack
option_chains_2x2                   = 0
option_chains_4x1                   = 0
option_ind_objects_mass             = 1
option_ind_objects_sfr              = 1
option_ind_objects_z                = 0 # 2 only, redshift dependence hack
option_samples_from_fitted_model    = 1
option_fitted_model_over_ms         = 1 # NA for mock
option_heatplots                    = 1 # legend labels not auto update

# =============================================================================
# inputs
# =============================================================================


massType = '_mStar'
#massType = '_mTot'

sfrType = '_delayed' # delayed means INSTANT
#sfrType = ''

massLim = '8p4'
zLim = '1p65'
wLim = '0p35'

#zLim = '2p5'
#wLim = '0p5'

chi2Lim = '9p5'
dimension = '2d'
comment = '201'
comment2    = '002'        # 001 is k==1, 002 is fitted k
nChains     = 4             # options are 1 or 4
minIter     = 5000            # options are 20, 200, 2000, or 
burn        = 4500

# =============================================================================
# load inputs
# =============================================================================

#sbf = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/linmix_AD_combined/linmix_npy_files_z{}_w{}_chi{}_dim{}_{}/{}_'.format(zLim, wLim, chi2Lim, dimension, comment, massLim)
    
sbf = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/linmix_AD_combined/linmix_npy_files_z{}_w{}_chi{}_dim{}{}{}_{}/{}_'.format(zLim, wLim, chi2Lim, dimension, massType, sfrType, comment, massLim)


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
    
GMMredshift        = np.load(sbf+'GMMredshift.npy')
GMMchi2     = np.load(sbf+'GMMchi2.npy')  
GMMmass     = np.load(sbf+'GMMmass.npy')  
GMMsfr      = np.load(sbf+'GMMsfr.npy')  

print(len(GMMmass))

# =============================================================================
# MAIN SEQUENCE
# =============================================================================

if option_main_sequence == 1:
 
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


#lm = np.load('/Users/lester/Documents/linmix_files/lm_chain_mass{}_z{}_w{}_chi{}_dim{}_com{}_{}x_{}_com{}.npy'.format(massLim, zLim, wLim, chi2Lim, dimension, comment, str(nChains), str(minIter), comment2))


lm = np.load('/Users/lester/Documents/linmix_files/lm_chain_mass{}_z{}_w{}_chi{}_dim{}{}{}_com{}_{}x_{}_com{}.npy'.format(massLim, zLim, wLim, chi2Lim, dimension, massType, sfrType, comment, str(nChains), str(minIter), comment2))

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
        plt.legend()
        plt.show()
        print(np.median(lm[param]), np.std(lm[param]))


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

if option_fitted_model_over_ms:

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
    
    z_low = float(zLim.replace('p','.')) - float(wLim.replace('p','.'))
    z_high = float(zLim.replace('p','.')) + float(wLim.replace('p','.'))
    
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), linewidth=lw, color='k', label=r'Fitted, {} $<$ z $<$ {}'.format(z_low, z_high))
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), linestyle='--', linewidth=lw, color='k', label=r'Fitted, Intrinsic Scatter')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), linestyle='--', linewidth=lw, color='k')
    
    plt.plot(s_x, s_y, color='w', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    plt.plot(s_x, s_y+s_y_intrinsic, color='w', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    plt.plot(s_x, s_y-s_y_intrinsic, color='w', linestyle='--', linewidth=lw)
    plt.plot(s_x, s_y+s_y_observed, color='w', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    plt.plot(s_x, s_y-s_y_observed, color='w', linestyle=':', linewidth=lw)
    
    plt.legend(loc='lower right')
    plt.show()    
    
    
    
    # DECENT PLOT - FIRST YEAR REPORT
    
    plt.figure(figsize=(10, 6))
    lw = 4
    #plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$M_\star$', size = 16)
    plt.ylabel(r'$\Psi$', size = 16)
    plt.hist2d(hp_x, hp_y, bins=50, range=((xlow, xhigh),(ylow, yhigh)), cmap=plt.cm.viridis)
    plt.colorbar(label='Count')
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax = plt.axes()
    ax.set_facecolor(rgba)
    
    z_low = float(zLim.replace('p','.')) - float(wLim.replace('p','.'))
    z_high = float(zLim.replace('p','.')) + float(wLim.replace('p','.'))
    
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), linewidth=lw, color='k', label=r'Fitted, {} $<$ z $<$ {}'.format(z_low, z_high))
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), linestyle='--', linewidth=lw, color='k', label=r'Fitted, Intrinsic Scatter')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), linestyle='--', linewidth=lw, color='k')
    
    plt.plot(s_x, s_y, color='w', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    plt.plot(s_x, s_y+s_y_intrinsic, color='w', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    plt.plot(s_x, s_y-s_y_intrinsic, color='w', linestyle='--', linewidth=lw)
    plt.plot(s_x, s_y+s_y_observed, color='w', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    plt.plot(s_x, s_y-s_y_observed, color='w', linestyle=':', linewidth=lw)
    
    plt.legend(loc='lower right')
    plt.tight_layout()
#    plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_results.png')
    plt.show()    

    # DECENT PLOT WITHOUT BOUNDS, also LOG
    
    plt.figure(figsize=(10, 6))
    lw = 4
    #plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
    plt.hist2d(hp_x, hp_y, bins=50, cmap=plt.cm.viridis, norm=mcolors.LogNorm())
#    plt.hist2d(hp_x, hp_y, bins=50, range=((xlow, xhigh),(ylow, yhigh)), cmap=plt.cm.viridis)
    plt.colorbar(label='Count')
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax = plt.axes()
#    ax.set_facecolor(rgba)
    
    z_low = float(zLim.replace('p','.')) - float(wLim.replace('p','.'))
    z_high = float(zLim.replace('p','.')) + float(wLim.replace('p','.'))
    
    plt.plot((xlow, xhigh), (alpha+beta*xlow, alpha+beta*xhigh), linewidth=lw, color='k', label=r'Fitted, {} $<$ z $<$ {}'.format(z_low, z_high))
    plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), linestyle='--', linewidth=lw, color='k', label=r'Fitted, Intrinsic Scatter')
    plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), linestyle='--', linewidth=lw, color='k')
    
    plt.plot(s_x, s_y, color='w', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    plt.plot(s_x, s_y+s_y_intrinsic, color='w', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    plt.plot(s_x, s_y-s_y_intrinsic, color='w', linestyle='--', linewidth=lw)
    plt.plot(s_x, s_y+s_y_observed, color='w', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    plt.plot(s_x, s_y-s_y_observed, color='w', linestyle=':', linewidth=lw)
    
    plt.legend(loc='lower right')
    plt.show()   

    # HEATPLOT POINTS    
    plt.scatter(hp_x, hp_y)
    plt.show()















