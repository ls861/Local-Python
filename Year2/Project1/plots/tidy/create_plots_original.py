=#!/usr/bin/env python2
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

load = True
load = False
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


# =============================================================================
# opening data
# =============================================================================
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3

if load:

    # converged chain
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_154_subset_002.p', 'rb') as f:
#    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_265_subset_001.p', 'rb') as f: # 9000-10000 
#        chain_MS = pickle.load(f, encoding='latin1') 
#    print(chain_burn.dtype.names)

    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_267_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_267 = pickle.load(f, encoding='latin1') #1.0
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
    with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_278_subset_001.p', 'rb') as f: # 9000-10000 
        chain_MS_278 = pickle.load(f, encoding='latin1') #k!=1




#    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p', 'rb') as f:
#        astrodeep_pickle = pickle.load(f, encoding='latin1') 
#    print(astrodeep_pickle.keys())

    # AD catalogue
    AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
    AD = np.load(AD_location)

    # from /Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/replicating_santini_with_santini_input_3dGMM.py
    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/data.p', 'rb') as f:
        data_dul = pickle.load(f, encoding='latin1') 
#    print(data_dul.keys())

    with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
        data = pickle.load(f, encoding='latin1') 
#    print(data.keys())

    bias_tests = ['110', '111', '112', '113', '114', '115', '116', '117', '118', '119']
    bias_tests = ['115', '117', '110', '119', '112', '111', '116', '113', '118', '114']
    chain_bias = []
    for test in bias_tests:
        with open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_truncated_{}.p'.format(test), 'rb') as f:
            chain_bias.append(pickle.load(f, encoding='latin1'))

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/mass_completeness_limits_0p90_new.fits'
    mass_completeness_limits_0p90_new = fits.open(fileName)
    #print(data_fits.info())
    #print(data_fits[1].header)

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.fits'
    s23 = fits.open(fileName)
    print(s23.info())
    print(s23[1].header)
    s23 = fits.open(fileName)[1].data # scenario 23


test = (s23['redshift_SANTINI']>4)
for i in range(len(s23['field_AD'][test])):
    print(s23['field_AD'][test][i])

for i in range(len(s23['field_AD'][test])):
    print(s23['id_BEAGLE'][test][i])
    

#%%    
# =============================================================================
# working out "errors" and medians on alpha, beta and ssfr
# =============================================================================

def ms_params(redshift, chain_MS):

    redshift_arr = np.repeat(np.array([redshift]).T, len(chain_MS), axis=1).T
    
    beta_a_arr = np.array([chain_MS['beta_a']]).T
    beta_b_arr = np.array([chain_MS['beta_b']]).T
    beta_arr = beta_a_arr + redshift_arr*beta_b_arr
    beta_16_arr = np.percentile(beta_arr, 16, axis=0)
#    beta_50_arr = np.median(beta_arr, axis=0)
    beta_84_arr = np.percentile(beta_arr, 84, axis=0)
    
    ssfr_a_arr = np.array([chain_MS['alphaN_a']]).T
    ssfr_b_arr = np.array([chain_MS['alphaN_b']]).T
    ssfr_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr)
    ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
#    ssfr_50_arr = np.median(ssfr_arr, axis=0)
    ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)
    
    alpha_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr) + normalisation - 9.0
    alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
#    alpha_50_arr = np.median(alpha_arr, axis=0)
    alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)
    
    # ARE THESE OK?
    beta_a = np.median(chain_MS['beta_a'])
    beta_b = np.median(chain_MS['beta_b'])
    ssfr_a = np.median(chain_MS['alphaN_a'])
    ssfr_b = np.median(chain_MS['alphaN_b'])
#    sig0          = np.median(chain_MS['sig0'])
#    pbad          = np.median(chain_MS['pbad'])
#    outlier_mean  = np.median(chain_MS['outlier_mean'])
#    outlier_sigma = np.median(chain_MS['outlier_sigma'])
    
    return ssfr_a, ssfr_b, beta_a, beta_b, ssfr_16_arr, ssfr_84_arr, beta_16_arr, beta_84_arr, alpha_16_arr, alpha_84_arr

normalisation = 9.7

z265 = np.linspace(1.25, 6.5, 1000)
z266 = np.linspace(1, 6.5, 1000)
z267 = np.linspace(1,2,1000)
z268 = np.linspace(2,3,1000)
z269 = np.linspace(3,4,1000)
z270 = np.linspace(4,5,1000)
z271 = np.linspace(1.25,6.5,1000)
z273 = np.linspace(1.25,2.0,1000)
z274 = np.linspace(1.25,6.5,1000)
z275 = np.linspace(1.25,2.0,1000)
z276 = np.linspace(1.25,2.0,1000)
z277 = np.linspace(1.25,2.0,1000)
z278 = np.linspace(1.25,2.0,1000)
z279 = np.linspace(1.25,2.0,1000)
z280 = np.linspace(1.25,2.0,1000)
z281 = np.linspace(1.25,6.5,1000)

ssfr_a_154, ssfr_b_154, beta_a_154, beta_b_154 = 0.056, 2.54, 0.687, -0.014 #154 original scenario 23 with Hogg etc
ssfr_a_217, ssfr_b_217, beta_a_217, beta_b_217 = 0.07, 2.351, 0.712, -0.075 #217 2x 3750 to 5000 without Hogg across all redshifts, -2<sfr<2
ssfr_a_218, ssfr_b_218, beta_a_218, beta_b_218 = 0.003, 6.831, 0.446, 0.234 #218 1st chain 4000burn
ssfr_a_219, ssfr_b_219, beta_a_219, beta_b_219 = 0.004, 6.44, 0.629, 0.139 #219 1st chain 4000burn median mass > 8.0
ssfr_a_220, ssfr_b_220, beta_a_220, beta_b_220 = 0.006, 6.435, 0.295, 0.589 #220 1st chain 4000burn median mass > 8.5
ssfr_a_221, ssfr_b_221, beta_a_221, beta_b_221 = 0.001, 8.583, -1.514, 2.011 #221 1st chain 4000burn median mass > 9.0
ssfr_a_232,ssfr_b_232,beta_a_232,beta_b_232 = 0.048,2.682,0.659,-0.003 #232 2x 3750 to 5000 ALL
ssfr_a_233,ssfr_b_233,beta_a_233,beta_b_233 = 0.137,1.854,0.766,-0.019 #233 2x 3750 to 5000 sfr cut, z>1
ssfr_a_234,ssfr_b_234,beta_a_234,beta_b_234 = 0.164,1.692,0.837,-0.053 #234 2x 3750 to 5000 z>1
ssfr_a_235,ssfr_b_235,beta_a_235,beta_b_235 = 0.315,1.209,0.975,-0.091 #235 2x 3750 to 5000 z>1.5
ssfr_a_259,ssfr_b_259,beta_a_259,beta_b_259 = 0.104,2.072,0.715,0.0 #259 2x 4500 to 6000 z>1.25, with Hogg, beta_b == 0
ssfr_a_263,ssfr_b_263,beta_a_263,beta_b_263 = 0.394,1.541,1.0,0.0 #263 2x 4500 to 6000 z>1.25, with Hogg, beta_b == 0, beta_a == 1, only nBurn covProposal
ssfr_a_264,ssfr_b_264,beta_a_264,beta_b_264 = 0.338,1.677,1.0,0.0 #264 2x 4500 to 6000 z>1.25, with Hogg, beta_b == 0, beta_a == 1, full covProposal
ssfr_a_266, ssfr_b_266, beta_a_266, beta_b_266 = 0.285, 1.802, 1.0, 0.0
ssfr_a_271, ssfr_b_271, beta_a_271, beta_b_271 = 0.104, 2.032, 0.712, 0.0
ssfr_a_273, ssfr_b_273, beta_a_273, beta_b_273 = 0.881, 0.0, 0.706, 0.0
ssfr_a_275, ssfr_b_275, beta_a_275, beta_b_275 = 1.003,0.0,0.725,0.0 #4x 9000-10000 mass 8.0, z 1.25-2.0
ssfr_a_276, ssfr_b_276, beta_a_276, beta_b_276 = 1.093,0.0,0.825,0.0 #4x 9000-10000 mass 8.5, z 1.25-2.0
ssfr_a_277, ssfr_b_277, beta_a_277, beta_b_277 = 7.925,0.0,1.83,0.0 #4x 9000-10000 mass 9.0, z 1.25-2.0
ssfr_a_278, ssfr_b_278, beta_a_278, beta_b_278 = 0.925,0.0,0.711,0.0 #4x 9000-10000 mass 8.0, z 1.25-2.0, k varied
ssfr_a_279, ssfr_b_279, beta_a_279, beta_b_279 = 1.205,0.0,0.844,0.0 #4x 9000-10000 mass 8.5, z 1.25-2.0, k varied
ssfr_a_280, ssfr_b_280, beta_a_280, beta_b_280 = 6.225,0.0,1.698,0.0 #4x 9000-10000 mass 9.0, z 1.25-2.0, k varied
ssfr_a_281, ssfr_b_281, beta_a_281, beta_b_281 = 1.625,0.0,0.836,0.0 #4x 9000-10000 mass 9.0, z 1.25-2.0

#ssfr_a_265, ssfr_b_265, beta_a_265, beta_b_265 = 0.375, 1.564, 1.0, 0.0
#ssfr_a_267, ssfr_b_267, beta_a_267, beta_b_267 = 0.723, 0.0, 0.714, 0.0
#ssfr_a_268, ssfr_b_268, beta_a_268, beta_b_268 = 1.316, 0.0, 0.797, 0.0
#ssfr_a_269, ssfr_b_269, beta_a_269, beta_b_269 = 2.217, 0.0, 0.818, 0.0
#ssfr_a_270, ssfr_b_270, beta_a_270, beta_b_270 = 4.035, 0.0, 0.472, 0.0
#ssfr_a_274, ssfr_b_274, beta_a_274, beta_b_274 = 0.132,1.975,0.773,0.0 

#ssfr_a_265, ssfr_b_265, beta_a_265, beta_b_265 = 0.375, 1.564, 1.0, 0.0
#ssfr_a_267, ssfr_b_267, beta_a_267, beta_b_267 = 0.723, 0.0, 0.714, 0.0
#ssfr_a_268, ssfr_b_268, beta_a_268, beta_b_268 = 1.316, 0.0, 0.797, 0.0
#ssfr_a_269, ssfr_b_269, beta_a_269, beta_b_269 = 2.217, 0.0, 0.818, 0.0
#ssfr_a_270, ssfr_b_270, beta_a_270, beta_b_270 = 4.035, 0.0, 0.472, 0.0
#ssfr_a_274, ssfr_b_274, beta_a_274, beta_b_274 = 0.132,1.975,0.773,0.0 

ssfr_a_265, ssfr_b_265, beta_a_265, beta_b_265, ssfr_16_arr_265, ssfr_84_arr_265, beta_16_arr_265, beta_84_arr_265, alpha_16_arr_265, alpha_84_arr_265 = ms_params(z265, chain_MS_265)#beta=1
ssfr_a_267, ssfr_b_267, beta_a_267, beta_b_267, ssfr_16_arr_267, ssfr_84_arr_267, beta_16_arr_267, beta_84_arr_267, alpha_16_arr_267, alpha_84_arr_267 = ms_params(z267, chain_MS_267)#1.0
ssfr_a_268, ssfr_b_268, beta_a_268, beta_b_268, ssfr_16_arr_268, ssfr_84_arr_268, beta_16_arr_268, beta_84_arr_268, alpha_16_arr_268, alpha_84_arr_268 = ms_params(z268, chain_MS_268)#2.0
ssfr_a_269, ssfr_b_269, beta_a_269, beta_b_269, ssfr_16_arr_269, ssfr_84_arr_269, beta_16_arr_269, beta_84_arr_269, alpha_16_arr_269, alpha_84_arr_269 = ms_params(z269, chain_MS_269)#3.0
ssfr_a_270, ssfr_b_270, beta_a_270, beta_b_270, ssfr_16_arr_270, ssfr_84_arr_270, beta_16_arr_270, beta_84_arr_270, alpha_16_arr_270, alpha_84_arr_270 = ms_params(z270, chain_MS_270)#4.0
ssfr_a_273, ssfr_b_273, beta_a_273, beta_b_273, ssfr_16_arr_273, ssfr_84_arr_273, beta_16_arr_273, beta_84_arr_273, alpha_16_arr_273, alpha_84_arr_273 = ms_params(z273, chain_MS_273)#1.25
ssfr_a_274, ssfr_b_274, beta_a_274, beta_b_274, ssfr_16_arr_274, ssfr_84_arr_274, beta_16_arr_274, beta_84_arr_274, alpha_16_arr_274, alpha_84_arr_274 = ms_params(z274, chain_MS_274)#beta=constant from z<4
ssfr_a_278, ssfr_b_278, beta_a_278, beta_b_278, ssfr_16_arr_278, ssfr_84_arr_278, beta_16_arr_278, beta_84_arr_278, alpha_16_arr_278, alpha_84_arr_278 = ms_params(z278, chain_MS_278)#k!=1


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
xhigh = 7.7

# =============================================================================
# beta
# =============================================================================

ax1 = fig.add_axes([0, 1, 0.5, 0.5]) #[left, bottom, width, height]

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

# redshift bins
#ax1.plot(z267, beta_a_267 + z267*beta_b_267, color='k') #1.0
ax1.plot(z268, beta_a_268 + z268*beta_b_268, color='k') #2.0
ax1.plot(z269, beta_a_269 + z269*beta_b_269, color='k') #3.0
ax1.plot(z270, beta_a_270 + z270*beta_b_270, color='k', label='z bins') #4.0
ax1.plot(z273, beta_a_273 + z273*beta_b_273, color='k') #1.25

# 1.25<z<6.5
#ax1.plot(z274, beta_a_154 + z274*beta_b_154, label='Original') # original
ax1.plot(z265, beta_a_265 + z265*beta_b_265, label='beta=1') # beta=1
#ax1.plot(z271, beta_a_271 + z271*beta_b_271, label='beta=constant') # beta=constant
ax1.plot(z274, beta_a_274 + z274*beta_b_274, label='z<4 for beta=constant') # beta=constant from z<4

# 1.25<z<2.0 for different lower masses
#ax1.plot(z275, beta_a_275 + z275*beta_b_275, label='8.0', color='red')
#ax1.plot(z276, beta_a_276 + z276*beta_b_276, label='8.5', color='lime')
#ax1.plot(z277, beta_a_277 + z277*beta_b_277, label='9.0', color='orange')
#ax1.plot(z278, beta_a_278 + z278*beta_b_278, label='8.0 with k', linestyle='dashed', color='red')
#ax1.plot(z279, beta_a_279 + z279*beta_b_279, label='8.5 with k', linestyle='dashed', color='lime')
#ax1.plot(z280, beta_a_280 + z280*beta_b_280, label='9.0 with k', linestyle='dashed', color='orange')

ax1.fill_between((0,0), (0,0), (0,0))
ax1.fill_between((0,0), (0,0), (0,0))
#ax1.fill_between(z267, beta_16_arr_267, beta_84_arr_267, alpha=0.3, color='k', zorder=0)
ax1.fill_between(z268, beta_16_arr_268, beta_84_arr_268, alpha=0.3, color='k', zorder=0)
ax1.fill_between(z269, beta_16_arr_269, beta_84_arr_269, alpha=0.3, color='k', zorder=0)
ax1.fill_between(z270, beta_16_arr_270, beta_84_arr_270, alpha=0.3, color='k', zorder=0)
ax1.fill_between(z273, beta_16_arr_273, beta_84_arr_273, alpha=0.3, color='k', zorder=0)
ax1.fill_between(z265, beta_16_arr_265, beta_84_arr_265, alpha=0.3, zorder=0)
ax1.fill_between(z274, beta_16_arr_274, beta_84_arr_274, alpha=0.3, zorder=0)


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

ax2 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]

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

# redshift bins
#ax2.plot(z267, np.log10(ssfr_a_267*(1.0+z267)**ssfr_b_267) + normalisation - 9.0, color='k') #1.0
ax2.plot(z268, np.log10(ssfr_a_268*(1.0+z268)**ssfr_b_268) + normalisation - 9.0, color='k') #2.0
ax2.plot(z269, np.log10(ssfr_a_269*(1.0+z269)**ssfr_b_269) + normalisation - 9.0, color='k') #3.0
ax2.plot(z270, np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0, color='k', label='z bins') #4.0
ax2.plot(z273, np.log10(ssfr_a_273*(1.0+z273)**ssfr_b_273) + normalisation - 9.0, color='k') #1.25

# 1.25<z<6.5
#ax2.plot(z274, np.log10(ssfr_a_154*(1.0+z274)**ssfr_b_154) + normalisation - 9.0, label='Original') # original
ax2.plot(z265, np.log10(ssfr_a_265*(1.0+z265)**ssfr_b_265) + normalisation - 9.0, label='beta=1') # beta=1
#ax2.plot(z271, np.log10(ssfr_a_271*(1.0+z271)**ssfr_b_271) + normalisation - 9.0, label='beta=constant') # beta=constant
ax2.plot(z274, np.log10(ssfr_a_274*(1.0+z274)**ssfr_b_274) + normalisation - 9.0, label='z<4 for beta=constant') # beta=constant from z<4

# 1.25<z<2.0 for different lower masses
#ax2.plot(z275, np.log10(ssfr_a_275*(1.0+z275)**ssfr_b_275) + normalisation - 9.0, label='8.0', color='red')
#ax2.plot(z276, np.log10(ssfr_a_276*(1.0+z276)**ssfr_b_276) + normalisation - 9.0, label='8.5', color='lime')
#ax2.plot(z277, np.log10(ssfr_a_277*(1.0+z277)**ssfr_b_277) + normalisation - 9.0, label='9.0', color='orange')
#ax2.plot(z278, np.log10(ssfr_a_278*(1.0+z278)**ssfr_b_278) + normalisation - 9.0, label='8.0 with k', linestyle='dashed', color='red')
#ax2.plot(z279, np.log10(ssfr_a_279*(1.0+z279)**ssfr_b_279) + normalisation - 9.0, label='8.5 with k', linestyle='dashed', color='lime')
#ax2.plot(z280, np.log10(ssfr_a_280*(1.0+z280)**ssfr_b_280) + normalisation - 9.0, label='9.0 with k', linestyle='dashed', color='orange')

ax2.fill_between((0,0), (0,0), (0,0))
ax2.fill_between((0,0), (0,0), (0,0))
#ax2.fill_between(z267, alpha_16_arr_267, alpha_84_arr_267, alpha=0.3, color='k', zorder=0)
ax2.fill_between(z268, alpha_16_arr_268, alpha_84_arr_268, alpha=0.3, color='k', zorder=0)
ax2.fill_between(z269, alpha_16_arr_269, alpha_84_arr_269, alpha=0.3, color='k', zorder=0)
ax2.fill_between(z270, alpha_16_arr_270, alpha_84_arr_270, alpha=0.3, color='k', zorder=0)
ax2.fill_between(z273, alpha_16_arr_273, alpha_84_arr_273, alpha=0.3, color='k', zorder=0)
ax2.fill_between(z265, alpha_16_arr_265, alpha_84_arr_265, alpha=0.3, zorder=0)
ax2.fill_between(z274, alpha_16_arr_274, alpha_84_arr_274, alpha=0.3, zorder=0)

ax2.set_xlim(xlow, xhigh)
ax2.set_xlabel('Redshift', labelpad=10)
ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax2.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(labelsize=fontsize_axes)

ax2.set_ylim(0.1, 1.9)
ax2.set_ylabel(string_normalisation, labelpad=10)
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax2.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(labelsize=fontsize_axes)

ax2.legend(bbox_to_anchor=(1, 0), loc=4, frameon=False, fontsize=fontsize_legend)

# =============================================================================
# ssfr
# =============================================================================

ax3 = fig.add_axes([0, 0, 0.5, 0.5]) #[left, bottom, width, height]
#ax3.plot(redshift, ssfr_50_arr, label='This work')
ax3.plot(z274, np.log10(ssfr_a_274*(1.0+z274)**ssfr_b_274)-9.0, label='This work')
ax3.fill_between(z274, ssfr_16_arr_274-9.0, ssfr_84_arr_274-9.0, alpha=0.3)

ax3.plot(z274, (np.log10((1+z274)**2.25)) + (np.log10(ssfr_a_274*(1+z274)**ssfr_b_274)) - np.log10((1+2.0)**2.25) -9.0 ) # normalise Dekel et al. 2009 to z=2

ax3.set_xlim(xlow, xhigh)
ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax3.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

ax3.set_ylim(-11.5, -6.5)
ax3.set_ylabel(string_ssfr, labelpad=10)
ax3.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax3.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(labelsize=fontsize_axes)

ax3.legend(bbox_to_anchor=(1, 0), loc=4, frameon=False, fontsize=fontsize_legend)

if save:
    plt.savefig('001_parameters_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#%%
# =============================================================================
# redshift vs redshift inc IRAC flags
# =============================================================================

idx1 = (data_dul['field_AD']%2.0==0.0) # clusters + parallels
idx2 = (data_dul['relflag_AD']==1.0) # relflag
idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

idx3_IRAC = np.logical_and(data_dul['CH1_BEAGLE_input']<-60.0, data_dul['CH2_BEAGLE_input']<-60.0)
idx3_IRAC = np.logical_and(idx3_IRAC, data_dul['redshift_BEAGLE']>4.0)
idx3_IRAC = ~idx3_IRAC

idx3_IRAC_removed = np.logical_and(data_dul['CH1_BEAGLE_input']<-60.0, data_dul['CH2_BEAGLE_input']<-60.0)
#idx3_IRAC_removed = np.logical_and(idx3_IRAC_removed, data_dul['Ks_BEAGLE_input']<-60.0)

MCLmassLow = np.empty(data_dul['redshift_BEAGLE'].shape)
for i in range(len(MCLmassLow)):
    if data_dul['redshift_BEAGLE'][i] <2.1789654:
        MCLmassLow[i] = 8.0
    elif data_dul['redshift_BEAGLE'][i] > 4.195:
        MCLmassLow[i] = 9.0
    else:
        MCLmassLow[i] = 6.91926521 + 0.49598529*data_dul['redshift_BEAGLE'][i]
idx4 = (data_dul['mass_BEAGLE_stellar'] + data_dul['mag_AD'] > MCLmassLow) 

idx5_1 = (abs(data_dul['redshift_BEAGLE']-((1.25+6.5)/2.0)) < (((1.25+6.5)/2.0) - 1.25)) # 0.5<redshift_BEAGLE<6.5 (only affects up to visual inspection of z=3.5)
idx5_2 = (abs(data_dul['redshift_BEAGLE']-data_dul['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
idx5_3 = np.logical_and((data_dul['redshift_BEAGLE'] < 3.5), (data_dul['redshift_AD'] < 3.5)) 
idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)   
idx5_zgt3p5 = np.full(len(data_dul['id_AD']), False)
for i in range(len(data_dul['id_AD'])):
    idx5_zgt3p5_temp = np.isin(data_dul['id_AD'][i],vis[:,1][(vis[:,0]==data_dul['field_AD'][i])&(vis[:,3]==1)])
    if idx5_zgt3p5_temp:
        idx5_zgt3p5[i] = True   
idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data_dul['redshift_BEAGLE'] < 6.5)
idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
   
TmassHigh = 9.244 + (0.753*np.minimum(data_dul['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data_dul['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
idx6 = (data_dul['mass_BEAGLE_stellar'] < (TmassHigh))

idx7 = (data_dul['new_min_chi2_BEAGLE']>0) & (data_dul['new_min_chi2_BEAGLE']<13.28) # chi-squared
idx8 = (data_dul['sfr_BEAGLE_instant']>-5.0) & (data_dul['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
idx9 = (data_dul['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103

# =============================================================================
# the plot (z vs z)
# =============================================================================

idx_contour = np.logical_and(idx2,idx3) #relflag + H
idx_contour = np.logical_and(idx_contour,idx7) #chi2
idx_IRAC = np.logical_and(idx_contour,idx3_IRAC_removed) # idx for no IRAC objects
idx_subset = np.logical_and(idx_contour,idx3_IRAC) # IRAC cut
idx_subset = np.logical_and(idx_subset,idx1) # clusters at the mo
idx_subset = np.logical_and(idx_subset,idx4) # lower mass cut
idx_subset = np.logical_and(idx_subset,idx6) # upper mass cut
idx_subset = np.logical_and(idx_subset,idx5_z) # redshift cuts + visual inspection
idx_subset = np.logical_and(idx_subset,idx8) # arbitrary sfr cut
idx_subset = np.logical_and(idx_subset,idx9) # GMM cut
print(sum(idx_subset))

###########
nbins_x = 30
nbins_y = 30
min_x = np.min(data_dul['redshift_AD'][idx_contour])
max_x = np.max(data_dul['redshift_AD'][idx_contour])
binsize_x = (max_x-min_x)/nbins_x
min_y = np.min(data_dul['redshift_BEAGLE'][idx_contour])
max_y = np.max(data_dul['redshift_BEAGLE'][idx_contour])
binsize_y = (min_y-max_y)/nbins_y
data_dul1 = plt.hist2d(data_dul['redshift_AD'][idx_contour], data_dul['redshift_BEAGLE'][idx_contour], bins=[40,40], range=[[min_x,max_x],[min_y,max_y]])
###########

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
ax1.contourf((data_dul1[1][1:]+data_dul1[1][:-1])/2., \
                (data_dul1[2][1:]+data_dul1[2][:-1])/2., np.log10(np.transpose(data_dul1[0])), \
                 cmap=cm.gist_yarg)
ax1.scatter(data_dul['redshift_AD'][idx_IRAC], data_dul['redshift_BEAGLE'][idx_IRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
ax1.scatter(data_dul['redshift_AD'][idx_subset], data_dul['redshift_BEAGLE'][idx_subset], s=3, alpha=1.0, color='red', label='Selected Subset')
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
      
#ax1.contourf((data_dul1[1][1:]+data_dul1[1][:-1])/2., \
#                (data_dul1[2][1:]+data_dul1[2][:-1])/2., np.log10(np.transpose(data_dul1[0])), \
#                 cmap=cm.gist_yarg)
#ax1.scatter(data_dul['redshift_AD'][idx_IRAC], data_dul['redshift_BEAGLE'][idx_IRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
#ax1.scatter(data_dul['redshift_AD'][idx_subset], data_dul['redshift_BEAGLE'][idx_subset], s=3, alpha=1.0, color='red', label='Selected Subset')


ax1.hist(((data_dul['redshift_BEAGLE']-data_dul['redshift_AD'])/(1.0+data_dul['redshift_AD']))[idx_subset], alpha=0.3, bins=40,log=True, label='median z')
#ax1.hist(((data_dul['redshift_BEAGLE_mean']-data_dul['redshift_AD'])/(1.0+data_dul['redshift_AD']))[idx_subset], alpha=0.3, bins=40, log=True, label='mean z')

#ax1.hist(data_dul['redshift_AD'][idx_subset], alpha=0.3, bins=30)
#ax1.hist(data_dul['redshift_BEAGLE'][idx_subset], alpha=0.3, bins=30)

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


ax1.scatter(data_dul['redshift_BEAGLE'][idx_subset], data_dul['mass_BEAGLE_stellar'][idx_subset] + data_dul['mag_AD'][idx_subset], s=3, alpha=1.0, color='blue', label='BEAGLE Fitted Mass', zorder=2)

ax1.scatter(data_dul['redshift_BEAGLE'][idx_subset], data_dul['mass_BEAGLE_stellar'][idx_subset], s=3, alpha=1.0, color='red', label='Magnification Corrected Mass', zorder=1)

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

ssfr_a, ssfr_b, beta_a, beta_b, sig0, pbad, outlier_mean, outlier_sigma = ssfr_a_274, ssfr_b_274, beta_a_274, beta_b_274, np.median(chain_MS_274['sig0']), np.median(chain_MS_274['pbad']), np.median(chain_MS_274['outlier_mean']), np.median(chain_MS_274['outlier_sigma'])

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

#ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')

x_tmp = np.array((0, 20))
ax1.plot(x_tmp, np.log10(ssfr_a*(1+2)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+3)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+4)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+5)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+6)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')

ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

idx_sort = np.argsort(data['redshift_BEAGLE'])
scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], data['sfr_BEAGLE_instant'][idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=1.25, vmax=6.5, alpha=1.0)

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

sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+normalisation-9.0-(normalisation*(beta_a+data['redshift_BEAGLE']*beta_b))

idx_sort = np.argsort(data['redshift_BEAGLE'])

scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=1.25, vmax=6.5, alpha=1.0)

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

sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+normalisation-9.0-(normalisation*(beta_a+data['redshift_BEAGLE']*beta_b))

log_p_xi_eta_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
p_bad = norm.pdf(data['sfr_BEAGLE_instant'][idx_sort], scale=outlier_sigma, loc=outlier_mean)

z_bad = pbad*p_bad
z_good = (1-pbad)*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

#ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
#ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')

x_tmp = np.array((0, 20))
ax1.plot(x_tmp, np.log10(ssfr_a*(1+2)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+3)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+4)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+5)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')
ax1.plot(x_tmp, np.log10(ssfr_a*(1+6)**ssfr_b) + normalisation - 9.0 + beta_a*(x_tmp-normalisation), color='k', alpha=0.5, linestyle='dashed')

ax1.plot((normalisation, normalisation), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], data['sfr_BEAGLE_instant'][idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5)

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

sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+normalisation-9.0-(normalisation*(beta_a+data['redshift_BEAGLE']*beta_b))

log_p_xi_eta_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
p_bad = norm.pdf(data['sfr_BEAGLE_instant'][idx_sort], scale=outlier_sigma, loc=outlier_mean)

z_bad = pbad*p_bad
z_good = (1-pbad)*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=np.log10((z_good/z_bad)[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-2.5, vmax=2.5)

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

idx_santini = ((data_dul['mass_BEAGLE_stellar']+data_dul['mag_AD'])[idx_subset] > 9.3) & (data_dul['redshift_BEAGLE'][idx_subset] > 1.0) & (data_dul['redshift_BEAGLE'][idx_subset] < 2.0)

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
#ax1.contourf((data_dul1[1][1:]+data_dul1[1][:-1])/2., \
#                (data_dul1[2][1:]+data_dul1[2][:-1])/2., np.log10(np.transpose(data_dul1[0])), \
#                 cmap=cm.gist_yarg)
#ax1.scatter(data_dul['redshift_AD'][idx_IRAC], data_dul['redshift_BEAGLE'][idx_IRAC], s=3, alpha=1.0, color='blue', label='No reliable IRAC data')
#ax1.scatter(data_dul['redshift_AD'][idx_subset], data_dul['redshift_BEAGLE'][idx_subset], s=3, alpha=1.0, color='red', label='Selected Subset')


ax1.hist((data_dul['mass_BEAGLE_stellar']-data_dul['mass_SANTINI'])[idx_subset][idx_santini], alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
ax1.hist((data_dul['sfr_BEAGLE_instant']-data_dul['sfr_SANTINI'])[idx_subset][idx_santini], alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])


#ax1.hist(((data_dul['mass_BEAGLE_stellar']-data_dul['mass_SANTINI'])/(1.0+data_dul['mass_SANTINI']))[idx_subset], alpha=0.3, bins=40,log=True, label='mass', range=[-2,2])
#ax1.hist(((data_dul['sfr_BEAGLE_instant']-data_dul['sfr_SANTINI'])/(1.0+data_dul['sfr_SANTINI']))[idx_subset], alpha=0.3, bins=40, log=True, label='sfr', range=[-2,2])

#ax1.hist(data_dul['redshift_AD'][idx_subset], alpha=0.3, bins=30)
#ax1.hist(data_dul['redshift_BEAGLE'][idx_subset], alpha=0.3, bins=30)

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
# Parameter vs redshift, in redshift bins
# =============================================================================

fig = plt.figure(figsize=(2*figuresize, 1*figuresize))
xlow = 0.3
xhigh = 5.7

# =============================================================================
# beta
# =============================================================================

ax1 = fig.add_axes([0, 1.0, 0.5, 0.5]) #[left, bottom, width, height]

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

ax2 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]

print(np.log10(ssfr_a_270*(1.0+z270)**ssfr_b_270) + normalisation - 9.0)

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

#print(np.median(chain_MS_267['pbad']), np.median(chain_MS_267['outlier_mean']), np.median(chain_MS_267['outlier_sigma']))
#print(np.median(chain_MS_268['pbad']), np.median(chain_MS_268['outlier_mean']), np.median(chain_MS_268['outlier_sigma']))
#print(np.median(chain_MS_269['pbad']), np.median(chain_MS_269['outlier_mean']), np.median(chain_MS_269['outlier_sigma']))
#print(np.median(chain_MS_270['pbad']), np.median(chain_MS_270['outlier_mean']), np.median(chain_MS_270['outlier_sigma']))
#
#print(len(chain_MS_267['pbad']))

ax3 = fig.add_axes([0, 0, 0.5, 0.5]) #[left, bottom, width, height]

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
ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax3.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

ax3.set_ylim(0.0, 0.65)
#ax3.set_ylim(0.3, 1.3)
ax3.set_ylabel(string_scatter, labelpad=10)
ax3.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax3.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(labelsize=fontsize_axes)

ax3.legend(bbox_to_anchor=(0.95, 0.05), loc=4, frameon=False, fontsize=fontsize_legend)

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








