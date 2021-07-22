#!/usr/bin/env python2
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
import pickle
from astropy.io import fits
import matplotlib.colors as mcolors
from scipy.stats import norm

# Collect all the font names available to matplotlib
#font_names = [f.name for f in fm.fontManager.ttflist]
#print(font_names)

# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
mpl.rc('image', cmap='jet')
cmap = mpl.cm.get_cmap('jet')
plt.rcParams['font.size'] = 24
plt.rcParams['axes.linewidth'] = 2
# plt.rcParams['text.usetex'] = True
fontsize_legend = 10
fontsize_axes = 24
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


#quick fix for windows
string_slope = r'slope'
string_normalisation = r'MS Normalisation'
string_scatter = r'scatter'
string_ssfr = r'ssfr'
string_mass = r'mass'
string_sfr = r'sfr'
string_deltaMS = r'deltaMS'
string_prob_ratio = r'prob ratio'
string_bias_test = r'bias test'
string_pbad = r'pbad'
string_outlier_mean = r'outlier mean'
string_outlier_sigma = r'outlier sigma'

# =============================================================================
# opening data
# =============================================================================
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3

def read_chain(zlow, zhigh, chain_MS, fit):
    
    z = np.linspace(zlow, zhigh, 1000)

    redshift_arr = np.repeat(np.array([z]).T, len(chain_MS['beta_a']), axis=1).T





    beta_a_arr = np.array([chain_MS['beta_a']]).T
    beta_b_arr = np.array([chain_MS['beta_b']]).T
#    print(np.shape(redshift_arr), np.shape(beta_b_arr))
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


# =============================================================================
# load data
# =============================================================================
    
if load:
    
    folder = "C:/Users/LSand/Documents/linmix_files/1000_sample_chains/"

    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_12_001.p", 'rb') as f:
        chain_MS_29_12_c1 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z2p0-3p0_4x10000_12_001.p", 'rb') as f:
        chain_MS_29_12_c2 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_12_001.p", 'rb') as f:
        chain_MS_29_12_c3 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z4p0-5p0_4x10000_12_001.p", 'rb') as f:
        chain_MS_29_12_c4 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z5p0-6p0_4x10000_12_001.p", 'rb') as f:
        chain_MS_29_12_c5 = pickle.load(f, encoding='latin1') 

    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z1p25-2p0_4x10000_14_001.p", 'rb') as f:
        chain_MS_29_14_c1 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z2p0-3p0_4x10000_14_001.p", 'rb') as f:
        chain_MS_29_14_c2 = pickle.load(f, encoding='latin1')  
    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z3p0-4p0_4x10000_14_001.p", 'rb') as f:
        chain_MS_29_14_c3 = pickle.load(f, encoding='latin1')   

    with open(folder+"PROCESSED_lm_chain_scenario_29_clusters_z3p0-4p0_4x50000_14_001.p", 'rb') as f:
        chain_MS_29_14_c3_50k = pickle.load(f, encoding='latin1')  

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
#    s29z1 = fits.open(fileName)
#    print(s29z1.info())
#    print(s29z1[1].header)
    s29z1 = fits.open(fileName)[1].data 

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z2p0-3p0.fits'
    s29z2 = fits.open(fileName)[1].data 

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z3p0-4p0.fits'
    s29z3 = fits.open(fileName)[1].data 

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z4p0-5p0.fits'
    s29z4 = fits.open(fileName)[1].data 

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z5p0-6p0.fits'
    s29z5 = fits.open(fileName)[1].data 
    '''
    rc_12_c1 = read_chain(1.25, 2.0, chain_MS_29_12_c1, 'const_beta_linear_alpha')
    rc_12_c2 = read_chain(2.0, 3.0, chain_MS_29_12_c2, 'const_beta_linear_alpha')
    rc_12_c3 = read_chain(3.0, 4.0, chain_MS_29_12_c3, 'const_beta_linear_alpha')
    rc_12_c4 = read_chain(4.0, 5.0, chain_MS_29_12_c4, 'const_beta_linear_alpha')
    rc_12_c5 = read_chain(5.0, 6.0, chain_MS_29_12_c5, 'const_beta_linear_alpha')

    rc_14_c1 = read_chain(1.25, 2.0, chain_MS_29_14_c1, 'const_beta_linear_alpha')
    rc_14_c2 = read_chain(2.0, 3.0, chain_MS_29_14_c2, 'const_beta_linear_alpha')
    rc_14_c3 = read_chain(3.0, 4.0, chain_MS_29_14_c3, 'const_beta_linear_alpha')

    rc_14_c3_50k = read_chain(3.0, 4.0, chain_MS_29_14_c3_50k, 'const_beta_linear_alpha')    
    

#%%
# =============================================================================
# redshift bin heatplots
# =============================================================================

s = [s29z1, s29z1, s29z1[(s29z1['mass_BEAGLE_stellar'] + s29z1['mag_AD']) > 8.5], s29z1[(s29z1['mass_BEAGLE_stellar'] + s29z1['mag_AD']) > 9.0], s29z2, s29z3, s29z4, s29z5, s29z3[(s29z3['sfr_BEAGLE_instant'] > -1.0) & (s29z3['sfr_BEAGLE_instant'] < 1.5)]]


rc = [rc_c1, rc_c1_lm8p0, rc_c1_lm8p5, rc_c1_lm9p0, rc_c2, rc_c3, rc_c4, rc_c5, rc_c3_nohogg]
z_lower = [1.25, 1.25, 1.25, 1.25, 2.0, 3.0, 4.0, 5.0, 3.0]
z_upper = [2.0, 2.0, 2.0, 2.0, 3.0, 4.0, 5.0, 6.0, 4.0]
chain = [chain_MS_29_c1, chain_MS_29_c1_lm8p0, chain_MS_29_c1_lm8p5, chain_MS_29_c1_lm9p0, chain_MS_29_c2, chain_MS_29_c3, chain_MS_29_c4, chain_MS_29_c5, chain_MS_29_c3_nohogg]

# need scenario 29 GMM values, and median values, both in ADx
# need to select just object with median within 1.25 - 2.0
# I think for paper I did a lower mass cut of 8? Do I want this still?

print(ADx.keys())

for f in range(len(s)):

    n_hp = 400 # number of samples to take from GMM in total
    
    x_hp = np.array([])
    y_hp = np.array([])
    z_hp = np.array([])
    
    for i in range(len(s[f]['mass_BEAGLE_stellar'])):
    
        for G in range(3):
            
            mean = np.array([s[f]['x_GMM_3d'][i,G],s[f]['y_GMM_3d'][i,G],s[f]['z_GMM_3d'][i,G]])
            cov = np.array([[np.power(s[f]['xsig_GMM_3d'][i,G],2), s[f]['xycov_GMM_3d'][i,G], s[f]['xzcov_GMM_3d'][i,G]],[s[f]['xycov_GMM_3d'][i,G], np.power(s[f]['ysig_GMM_3d'][i,G],2), s[f]['yzcov_GMM_3d'][i,G]],[s[f]['xzcov_GMM_3d'][i,G], s[f]['yzcov_GMM_3d'][i,G], np.power(s[f]['zsig_GMM_3d'][i,G],2)]])
            
            xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s[f]['amp_GMM_3d'][i,G]))
    
            x_hp = np.concatenate((x_hp,xyz[:,0]))
            y_hp = np.concatenate((y_hp,xyz[:,1]))
            z_hp = np.concatenate((z_hp,xyz[:,2]))
            
    # only keep GMM samples within the redshift bin
    x_hp = x_hp[abs(z_hp - ((z_upper[f]+z_lower[f])/2.0)) < ((z_upper[f]-z_lower[f])/2.0)]
    y_hp = y_hp[abs(z_hp - ((z_upper[f]+z_lower[f])/2.0)) < ((z_upper[f]-z_lower[f])/2.0)]
    z_hp = z_hp[abs(z_hp - ((z_upper[f]+z_lower[f])/2.0)) < ((z_upper[f]-z_lower[f])/2.0)]
    
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
    #ax1.plot((8.5,8.5), (ylow, yhigh), color='w', linestyle='dashed', linewidth=2)
    #ax1.plot((10.0,10.0), (ylow, yhigh), color='w', linestyle='dashed', linewidth=2)
    
    x_tmp = np.array([xlow, xhigh])
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, {} $<$ z $<$ {}'.format(str(z_lower[f]), str(z_upper[f])))
    #ax1.plot(x_tmp, (x_tmp-9.7)*rc_c1['beta_50_arr'][0] + rc_c1['alpha_50_arr'][0], color='w')
    
    xi_min = 8.5
    xi_max = 10.0    
    sig = rc[f]['sig0_50_arr'][0] * ( ((1.0-rc[f]['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)
    
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
    
    #ax1.plot(s_x, s_y, color='k', linewidth=lw, label=r'Santini+17, 1.3 $<$ z $<$ 2.0')
    #ax1.plot(s_x, s_y+s_y_intrinsic, color='k', linestyle='--', linewidth=lw, label=r'Santini+17, Intrinsic Scatter')
    #ax1.plot(s_x, s_y-s_y_intrinsic, color='k', linestyle='--', linewidth=lw)
    #ax1.plot(s_x, s_y+s_y_observed, color='k', linestyle=':', linewidth=lw, label=r'Santini+17, Observed Scatter')
    #ax1.plot(s_x, s_y-s_y_observed, color='k', linestyle=':', linewidth=lw)
    
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
    

    # =============================================================================
    # redshift bin MS colour coded by MS probability
    # =============================================================================
    
    medians = get_medians(chain[f])
    
    #sfr_surface_real = ((medians['beta_a']+s['redshift_BEAGLE']*medians['beta_b'])*s['mass_BEAGLE_stellar'])+(np.log10(medians['alphaN_a']*((1+s['redshift_BEAGLE'])**medians['alphaN_b'])))+normalisation-9.0-(normalisation*(medians['beta_a']+s['redshift_BEAGLE']*medians['beta_b']))
    
    sfr_surface_real = ((medians['beta_a']+s[f]['redshift_BEAGLE']*medians['beta_b'])*(s[f]['mass_BEAGLE_stellar'] - 9.7))  + \
                        (medians['alphaN_a']+s[f]['redshift_BEAGLE']*medians['alphaN_b'])
    
#    log_p_xi_eta_theta = norm.logpdf(s[f]['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
    log_p_eta_xi_theta = norm.logpdf(s[f]['sfr_BEAGLE_instant'], scale=medians['sig0'], loc=sfr_surface_real)
    p_bad = norm.pdf(s[f]['sfr_BEAGLE_instant'], scale=medians['outlier_sigma'], loc=medians['outlier_mean'])
    
    z_bad = medians['pbad']*p_bad
    z_good = (1.0-medians['pbad'])*np.exp(log_p_eta_xi_theta)
    
    z_bad = np.where(z_bad==0, 1e-9, z_bad) # necessary if NO HOGG, ie pbad == 0
    
    idx_sort = np.argsort(z_good/z_bad)
    
    fig = plt.figure(figsize=(figuresize, figuresize))
    ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
    
    
    ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')
    
    scatter = ax1.scatter(s[f]['mass_BEAGLE_stellar'][idx_sort], s[f]['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5, s=100)
    
    xlow = 6.5
    xhigh = 10.5
    ylow = -2.5
    yhigh = 2.5
    
    ximin = 8.5
    ximax = 10.0
    
    lw = 3
    
    ax1.plot((9.7,9.7), (ylow, yhigh), color='gray', linestyle='dashed', linewidth=2)
    
    x_tmp = np.array([xlow, xhigh])
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0], color='r', linewidth=lw, label=r'Our Work, {} $<$ z $<$ {}'.format(str(z_lower[f]), str(z_upper[f])))
     
    sig = rc[f]['sig0_50_arr'][0] * ( ((1.0-rc[f]['k_50_arr'][0])*(x_tmp-xi_max)/(xi_max-xi_min)) + 1.0 )
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0] + sig, color='r', linestyle='dashed', linewidth=lw, label=r'Our Work, Intrinsic Scatter')
    ax1.plot(x_tmp, (x_tmp-9.7)*rc[f]['beta_50_arr'][0] + rc[f]['alpha_50_arr'][0] - sig, color='r', linestyle='dashed', linewidth=lw)
    
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
    



#%%

# =============================================================================
# Parameter vs redshift
# =============================================================================

rcs = [rc_12_c1, rc_12_c2, rc_12_c3, rc_12_c4, rc_12_c5, rc_14_c1, rc_14_c2, rc_14_c3, rc_14_c3_50k]

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = -0.3
xhigh = 1.3 + len(rcs) + 2

param = ['beta', 'alpha', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']

ax1 = fig.add_axes([0, 2.5, 0.5, 0.5]) #[left, bottom, width, height]
ax2 = fig.add_axes([0, 2.0, 0.5, 0.5]) #[left, bottom, width, height]
ax3 = fig.add_axes([0, 1.5, 0.5, 0.5]) #[left, bottom, width, height]
ax4 = fig.add_axes([0.6, 2.5, 0.5, 0.5]) #[left, bottom, width, height]
ax5 = fig.add_axes([0.6, 2.0, 0.5, 0.5]) #[left, bottom, width, height]
ax6 = fig.add_axes([0.6, 1.5, 0.5, 0.5]) #[left, bottom, width, height]

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


for rc in rcs:

    for a, ax in enumerate(axes):
        ax.plot(rc['z']-min(rc['z'])+count, rc['{}_50_arr'.format(param[a])])
        ax.fill_between(rc['z']-min(rc['z'])+count, rc['{}_16_arr'.format(param[a])], rc['{}_84_arr'.format(param[a])], alpha=0.3, zorder=0)
    count+=1

for a, ax in enumerate(axes):
    ax.plot(0, 0, color='#1f77b4', label='12 z1')
    ax.plot(0, 0, color='#ff7f0e', label='12 z2')
    ax.plot(0, 0, color='#2ca02c', label='12 z3')
    ax.plot(0, 0, color='#d62728', label='12 z4')
    ax.plot(0, 0, color='#9467bd', label='12 z5')
    ax.plot(0, 0, color='#8c564b', label='14 z1')
    ax.plot(0, 0, color='#e377c2', label='14 z2')
    ax.plot(0, 0, color='gray', label='14 z3')
    ax.plot(0, 0, color='yellow', label='14 z3 50k')
    ax.legend(loc='upper right', frameon=False, fontsize=fontsize_legend)

# layout
ylim_low = np.array([0.2, 0.0, 0.0, 0.0, -6.0, 0.0])
ylim_high = np.array([1.4, 2.5, 1.0, 1.0, 6.0, 6.0])
string = [string_slope, string_normalisation, string_scatter, string_pbad, string_outlier_mean, string_outlier_sigma]
# ticker_maj = (ylim_high - ylim_low) / 3.0
ticker_maj = np.array([0.4, 0.5, 0.2, 0.2, 3, 2]) 
ticker_min = ticker_maj / 2.0
for a, ax in enumerate(axes):
    ax.set_ylim(ylim_low[a], ylim_high[a])
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












