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
fontsize_legend = 18
fontsize_axes = 24
figuresize = 7

normalisation = 9.7

#load = True
load = False
save = False
 
# =============================================================================
# useful axis strings
# =============================================================================

string_slope = r'$\mathrm{MS\,Slope,}\,\beta$'
string_normalisation = r'MS Normalisation, $\alpha$'
string_scatter = r'$\mathrm{Scatter \, around \, the \, MS}$'
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


if load:

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
    s29z1 = fits.open(fileName)
#    print(s29z1.info())
#    print(s29z1[1].header)
    s29z1 = fits.open(fileName)[1].data # scenario 23

    
    rc_c1 = read_chain(1.25, 2.0, chain_MS_29_c1, 'const_beta_ssfr_alpha')
    rc_c2 = read_chain(2.0, 3.0, chain_MS_29_c2, 'const_beta_ssfr_alpha')
    rc_c3 = read_chain(3.0, 4.0, chain_MS_29_c3, 'const_beta_ssfr_alpha')
    rc_c4 = read_chain(4.0, 5.0, chain_MS_29_c4, 'const_beta_ssfr_alpha')
    rc_c5 = read_chain(5.0, 6.0, chain_MS_29_c5, 'const_beta_ssfr_alpha')
    rc_c = read_chain(1.25, 6.0, chain_MS_29_c, 'const_beta_ssfr_alpha')
    rc_c1k = read_chain(1.25, 2.0, chain_MS_29_c1k, 'const_beta_ssfr_alpha')




#%%
# =============================================================================
# mass dependent scatter heatplot 1.25<z<2.0
# =============================================================================


# need scenario 29 GMM values, and median values, both in ADx
# need to select just object with median within 1.25 - 2.0
# I think for paper I did a lower mass cut of 8? Do I want this still?
    
print(ADx.keys())

n_hp = 100 # number of samples to take from GMM in total

x_hp = np.array([])
y_hp = np.array([])
z_hp = np.array([])

for i in range(len(s29z1['mass_BEAGLE_stellar'])):
    
    for G in range(3):
        
        mean = np.array([s29z1['x_GMM_3d'][i,G],s29z1['y_GMM_3d'][i,G],s29z1['z_GMM_3d'][i,G]])
        cov = np.array([[np.power(s29z1['xsig_GMM_3d'][i,G],2), s29z1['xycov_GMM_3d'][i,G], s29z1['xzcov_GMM_3d'][i,G]],[s29z1['xycov_GMM_3d'][i,G], np.power(s29z1['ysig_GMM_3d'][i,G],2), s29z1['yzcov_GMM_3d'][i,G]],[s29z1['xzcov_GMM_3d'][i,G], s29z1['yzcov_GMM_3d'][i,G], np.power(s29z1['zsig_GMM_3d'][i,G],2)]])

        xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s29z1['amp_GMM_3d'][i,G]))

        x_hp = np.concatenate((x_hp,xyz[:,0]))
        y_hp = np.concatenate((y_hp,xyz[:,1]))
        z_hp = np.concatenate((z_hp,xyz[:,2]))

# only keep GMM samples within the redshift bin
x_hp = x_hp[abs(z_hp - 1.625) < 0.375]
y_hp = y_hp[abs(z_hp - 1.625) < 0.375]
z_hp = z_hp[abs(z_hp - 1.625) < 0.375]

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
#    ax1.set_title('{} - {}'.format(z_med_hp_low, z_med_hp_high))
lw = 3

xlow = 5.5
xhigh = 11.5
ylow = -3.5
yhigh = 3.5

xlow = 6.5
xhigh = 10.5
ylow = -2.5
yhigh = 2.5

ximin = 8.5
ximax = 10.0

h = ax1.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)

#ax1.plot((8.0,8.0), (ylow, yhigh), color='w', linestyle='dashed', linewidth=2)
ax1.plot((8.5,8.5), (ylow, yhigh), color='w', linestyle='dashed', linewidth=2)
ax1.plot((10.0,10.0), (ylow, yhigh), color='w', linestyle='dashed', linewidth=2)


x_tmp = np.array([xlow, xhigh])
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, 1.25 $<$ z $<$ 2.0')
#ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1['beta_50_arr'][0] + rc_c1['alpha_50_arr'][0], color='w')

xi_min = 8.5
xi_max = 10.0    
sig = rc_c1k['sig0_50_arr'][0] * ( ((1.0-rc_c1k['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1k['beta_50_arr'][0] + rc_c1k['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)

'''
print(rc_c1['alpha_50_arr'][0])
print(rc_c1k['alpha_50_arr'][0])
'''

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
    
'''
print('alpha', s1_alpha*(0 - 9.7) + s1_beta)
print('beta', (s1_alpha*(1 - 9.7) + s1_beta) - (s1_alpha*(0 - 9.7) + s1_beta) )
'''

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

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(h[3], cax=cbaxes)
cb.set_ticks([5, 10, 50, 100])
cb.set_ticklabels([5, 10, 50, 100])
#cbaxes.set_ylabel(r'TEST', rotation=270, labelpad=30)
#cbaxes.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=1))
cbaxes.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
cbaxes.yaxis.set_tick_params(which='minor', size=0)
cbaxes.yaxis.set_tick_params(labelsize=fontsize_axes)

cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax1.set_facecolor(rgba)

ax1.legend(bbox_to_anchor=(1, 0), loc='lower right', frameon=True, fontsize=0.8*fontsize_legend)

if save:
    plt.savefig('KICC_1.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()




#%%
# =============================================================================
# sigma at mass 8.5 and 10.0 histograms
# =============================================================================

fig = plt.figure(figsize=(0.8*figuresize, 0.8*figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]

ax1.hist(chain_MS_29_c1k['sig0']*chain_MS_29_c1k['k'], alpha=0.3, bins=30, density=True, range=[0.1, 0.5], label=r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 8.5$')
ax1.hist(chain_MS_29_c1k['sig0'], alpha=0.3, bins=30, density=True, range=[0.1, 0.5], label=r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$')

'''
$\mathrm{Stellar} \,  \mathrm{Mass} = 10^{8.5}\, \mathrm{M_{\odot}}$
$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 8.5$

$\mathrm{Stellar} \,  \mathrm{Mass} = 10^{10.0}\, \mathrm{M_{\odot}}$
$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$
'''

#ax1.hist(((ADx_subset['mass_BEAGLE_stellar']-ADx_subset['mass_SANTINI'])/(1.0+ADx_subset['mass_SANTINI'])), alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
#ax1.hist(((ADx_subset['sfr_BEAGLE_instant']-ADx_subset['sfr_SANTINI'])/(1.0+ADx_subset['sfr_SANTINI'])), alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])

#ax1.hist(ADx_subset['redshift_AD'], alpha=0.3, bins=30)
#ax1.hist(ADx_subset['redshift_BEAGLE'], alpha=0.3, bins=30)

#ax1.set_xlim(-50, 50)
ax1.set_xlabel(r'$\mathrm{Intrinsic} \,  \mathrm{Scatter}$', labelpad=10)
#ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.set_ylim(0.0, 10.0)
#ax1.set_ylabel(r'Count', labelpad=10)
#ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=0, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=0, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=0*fontsize_axes)


ax1.legend(loc='lower right', frameon=True, fontsize=0.8*fontsize_legend)

if save:
    plt.savefig('KICC_2.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#%%
# =============================================================================
# ssfr and scatter vs redshift
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1.5*figuresize))
xlow = -0.3
xhigh = 7.3

param = ['ssfr', 'sig0']

ax1 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 0.0, 0.5, 0.5]) #[left, bottom, width, height]

axes = [ax1,ax2]

for ax in [ax1]:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=0, width=2, direction='out', bottom='on', top='on')
    ax.xaxis.set_tick_params(which='minor', size=0, width=2, direction='out', bottom='on', top='on')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=0, width=2, direction='out')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)

for ax in [ax2]:
    ax.set_xlim(xlow, xhigh)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.xaxis.set_tick_params(which='minor', size=0, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=10, width=2, direction='out')
    ax.yaxis.set_tick_params(which='minor', size=0, width=2, direction='out')
    ax.yaxis.set_tick_params(labelsize=fontsize_axes)
    
    
# redshift bins
for rc in [rc_c1, rc_c2, rc_c3, rc_c4, rc_c5]:
    for a, ax in enumerate([ax2]):
        a += 1
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='red', linewidth=lw)
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='red', zorder=0)

for rc in [rc_c]:
    for a, ax in enumerate([ax1]):
        ax.plot(rc['z'], rc['{}_50_arr'.format(param[a])], color='red', linewidth=lw)
        ax.fill_between(rc['z'], rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, color='red', zorder=0)
        

     
# layout
ylim_low = [-10.5, -0.1]
ylim_high = [-7.3, 0.79]
string = [string_ssfr, string_scatter]
ticker_maj = [1.0, 0.2]
ticker_min = [0.1, 0.1]
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high [a])
    ax.set_ylabel(string[a], labelpad=10)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ticker_maj[a]))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ticker_min[a]))

ax2.set_xlabel('Redshift', labelpad=10)
ax2.xaxis.set_tick_params(labelsize=fontsize_axes)

#ax2.plot(0, 0, color='red', label='Our Work', linewidth=lw)
ax2.plot(0, 0, color='red', label='Our Work, Redshift Bins', linewidth=lw)

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

ax1.legend(loc='lower right', frameon=False, fontsize=0.7*fontsize_legend, ncol=2)

plt.text(0.0, 1.6, r'$\log(\mathrm{Stellar} \, \mathrm{Mass} \, / \, \mathrm{M_{\odot}}) \sim 9.7$')

# santini scatter
ax2.plot((1.3,2.0), (0.0,0.0), color='k', linewidth=lw, label=r'Santini+17, Intrinsic Scatter,  $\log(\mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$')
ax2.plot((2.0,3.0), (0.0,0.0), color='k', linewidth=lw)
ax2.plot((3.0,4.0), (0.33,0.33), color='k', linewidth=lw)
ax2.plot((4.0,5.0), (0.04,0.04), color='k', linewidth=lw)

ax2.plot((1.3,2.0), (0.26,0.26), color='k', linewidth=lw, linestyle='--', label=r'Santini+17, Observed Scatter, $\log(\mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 10.0$')
ax2.plot((2.0,3.0), (0.22,0.22), color='k', linewidth=lw, linestyle='--')
ax2.plot((3.0,4.0), (0.45,0.45), color='k', linewidth=lw, linestyle='--')
ax2.plot((4.0,5.0), (0.46,0.46), color='k', linewidth=lw, linestyle='--')

# santini scatter
ax2.plot((1.3,2.0), (0.36,0.36), color='k', linewidth=lw, linestyle='dotted', label=r'Santini+17, Intrinsic Scatter, $\log(\mathrm{Mass} \, / \, \mathrm{M_{\odot}}) = 8.4$')
ax2.plot((2.0,3.0), (0.27,0.27), color='k', linewidth=lw, linestyle='dotted')
ax2.plot((3.0,4.0), (0.0,0.0), color='k', linewidth=lw, linestyle='dotted')
ax2.plot((4.0,5.0), (0.34,0.34), color='k', linewidth=lw, linestyle='dotted')



ax2.legend(loc='upper left', frameon=False, fontsize=0.7*fontsize_legend)


if save:
    plt.savefig('KICC_3.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

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

medians = get_medians(chain_MS_29_c1k)

#from astropy.cosmology import FlatLambdaCDM
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
#
#z_med_spe = np.array([0.0, 1.0, 6.0, 9999.0])
#t_med_spe = cosmo.age(z_med_spe).value
#print(t_med_spe)

#%%
# =============================================================================
# MS colour coded by MS probability
# =============================================================================

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

ax1.legend(loc='lower right', frameon=True, fontsize=0.8*fontsize_legend, framealpha=1.0)


cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)

cb.set_ticks(np.linspace(-2, 2, 5))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\mathrm{MS}\longleftarrow \,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\, \longrightarrow \mathrm{Outliers}$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=0, width=2, direction='in', labelsize=0*fontsize_axes)

if save:
    plt.savefig('KICC_4.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()





#%%
# =============================================================================
# combining 2 and 4
# =============================================================================





fig = plt.figure(figsize=(0.8*figuresize, 0.8*figuresize))


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

ax1.legend(loc='lower right', frameon=True, fontsize=0.8*fontsize_legend, framealpha=1.0)


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


ax2.legend(loc='lower right', frameon=True, fontsize=0.8*fontsize_legend)




#save = True
if save:
    plt.savefig('KICC_5.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()
















