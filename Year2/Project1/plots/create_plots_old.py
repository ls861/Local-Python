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
fontsize_legend = 20
fontsize_axes = 20
figuresize = 7

# =============================================================================
# slope, normalisation and ssfr
# =============================================================================

ssfr_a = 0.056
ssfr_b = 2.54
beta_a = 0.687
beta_b = -0.014

redshift = np.linspace(0.3, 6.5, 1000)
beta = beta_a + redshift*beta_b
alpha = np.log10(ssfr_a*(1.0+redshift)**ssfr_b) + 9.7 - 9.0
ssfr = np.log10(ssfr_a*(1.0+redshift)**ssfr_b)


# =============================================================================
# working out "errors"
# =============================================================================

'''
with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_154_subset_002.p', 'rb') as f:
    chain_burn = pickle.load(f, encoding='latin1') 
'''
#print(chain_burn.dtype.names)
    
redshift_arr = np.repeat(np.array([redshift]).T, len(chain_burn), axis=1).T
beta_a_arr = np.array([chain_burn['beta_a']]).T
beta_b_arr = np.array([chain_burn['beta_b']]).T
beta_arr = beta_a_arr + redshift_arr*beta_b_arr
beta_16_arr = np.percentile(beta_arr, 16, axis=0)
beta_84_arr = np.percentile(beta_arr, 84, axis=0)
#beta_1_arr = np.percentile(beta_arr, 1, axis=0)
#beta_99_arr = np.percentile(beta_arr, 99, axis=0)

ssfr_a_arr = np.array([chain_burn['alphaN_a']]).T
ssfr_b_arr = np.array([chain_burn['alphaN_b']]).T
alpha_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr) + 9.7 - 9.0
alpha_16_arr = np.percentile(alpha_arr, 16, axis=0)
alpha_84_arr = np.percentile(alpha_arr, 84, axis=0)

ssfr_arr = np.log10(ssfr_a_arr*(1.0+redshift_arr)**ssfr_b_arr)
ssfr_16_arr = np.percentile(ssfr_arr, 16, axis=0)
ssfr_84_arr = np.percentile(ssfr_arr, 84, axis=0)

normalisation = 9.7


# =============================================================================
# 154 and 217 redshift bin test
# =============================================================================

test_z = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])

#154 values at 7.7 and 9.7 +- 0.1 np.mean (z from 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
test_154_7p7 = np.array([-4.31,-0.770,-0.227,0.371,1.02,1.08])
test_154_8p7 = np.array([-3.21,-0.663,-0.0991,0.380,1.30,1.37])
test_154_9p7 = np.array([-1.56,0.684,0.457,0.801,1.753,1.81])
test_154_beta = (test_154_9p7 - test_154_7p7) / 2.0
test_154_beta_1 = test_154_9p7 - test_154_8p7
test_154_beta_2 = test_154_8p7 - test_154_7p7

#217 values at 7.7 and 9.7 +- 0.1 np.mean  (z from 0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
test_217_7p7 = np.array([-2.42,-0.559,-0.199,0.330,0.993,1.12])
test_217_8p7 = np.array([-1.52,-0.277,-0.120,0.441,1.18,1.24])
test_217_9p7 = np.array([-0.131,0.491,0.754,0.811,1.51,1.43])
test_217_beta = (test_217_9p7 - test_217_7p7) / 2.0
test_217_beta_1 = test_217_9p7 - test_217_8p7
test_217_beta_2 = test_217_8p7 - test_217_7p7

print(test_154_beta)
print(test_217_beta)
print(test_154_9p7)
print(test_217_9p7)

#217 
ssfr_a_217 = 0.07
ssfr_b_217 = 2.351
beta_a_217 = 0.712
beta_b_217 = -0.075

redshift_217 = np.linspace(0.3, 6.5, 1000)
beta_217 = beta_a_217 + redshift_217*beta_b_217
alpha_217 = np.log10(ssfr_a_217*(1.0+redshift_217)**ssfr_b_217) + 9.7 - 9.0
ssfr_217 = np.log10(ssfr_a_217*(1.0+redshift_217)**ssfr_b_217)

# HOGG outlier sigma limited to > 2

#228 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.048,2.724,0.681,-0.033
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_228 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_228 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_228 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#229 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.221,1.485,0.894,-0.098
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_229 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_229 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_229 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#230 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.193,1.612,0.873,-0.071
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_230 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_230 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_230 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#231 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.372,1.178,0.936,-0.087
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_231 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_231 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_231 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

# HOGG BACK TO NORMAL

#232 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.048,2.682,0.659,-0.003
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_232 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_232 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_232 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#233 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.137,1.854,0.766,-0.019
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_233 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_233 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_233 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#234 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.164,1.692,0.837,-0.053
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_234 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_234 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_234 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#235 2x 3750 to 5000
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.315,1.209,0.975,-0.091
redshift_tmp = np.linspace(0.3, 6.5, 1000)
beta_235 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_235 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_235 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)



# REDSHIFT BINS

#236 2250 - 3000 x2
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.021,3.698,0.67,-0.086
redshift_236 = np.linspace(0.0, 1.0, 1000)
beta_236 = beta_a_tmp + redshift_236*beta_b_tmp
alpha_236 = np.log10(ssfr_a_tmp*(1.0+redshift_236)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_236 = np.log10(ssfr_a_tmp*(1.0+redshift_236)**ssfr_b_tmp)

#237 2250 - 3000 x2
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.207,1.257,0.942,-0.186
redshift_237 = np.linspace(1.0, 1.5, 1000)
beta_237 = beta_a_tmp + redshift_237*beta_b_tmp
alpha_237 = np.log10(ssfr_a_tmp*(1.0+redshift_237)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_237 = np.log10(ssfr_a_tmp*(1.0+redshift_237)**ssfr_b_tmp)

#238 2250 - 3000 c0
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.01,5.876,0.888,0.041
redshift_238 = np.linspace(1.5, 2.0, 1000)
beta_238 = beta_a_tmp + redshift_238*beta_b_tmp
alpha_238 = np.log10(ssfr_a_tmp*(1.0+redshift_238)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_238 = np.log10(ssfr_a_tmp*(1.0+redshift_238)**ssfr_b_tmp)

#239 2250 - 3000 x2
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.16,1.821,0.198,0.276
redshift_239 = np.linspace(2.0, 2.5, 1000)
beta_239 = beta_a_tmp + redshift_239*beta_b_tmp
alpha_239 = np.log10(ssfr_a_tmp*(1.0+redshift_239)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_239 = np.log10(ssfr_a_tmp*(1.0+redshift_239)**ssfr_b_tmp)

#240 2250 - 3000 c0
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.475,0.452,0.331,0.118
redshift_240 = np.linspace(2.5, 3.0, 1000)
beta_240 = beta_a_tmp + redshift_240*beta_b_tmp
alpha_240 = np.log10(ssfr_a_tmp*(1.0+redshift_240)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_240 = np.log10(ssfr_a_tmp*(1.0+redshift_240)**ssfr_b_tmp)

#241a 2250 - 3000 c0
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.001,5.125,0.124,0.179
redshift_241a = np.linspace(3.0, 6.5, 1000)
beta_241a = beta_a_tmp + redshift_241a*beta_b_tmp
alpha_241a = np.log10(ssfr_a_tmp*(1.0+redshift_241a)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_241a = np.log10(ssfr_a_tmp*(1.0+redshift_241a)**ssfr_b_tmp)

#241b 2250 - 3000 c1
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.114,2.002,1.109,-0.104
redshift_241b = np.linspace(3.0, 6.5, 1000)
beta_241b = beta_a_tmp + redshift_241b*beta_b_tmp
alpha_241b = np.log10(ssfr_a_tmp*(1.0+redshift_241b)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_241b = np.log10(ssfr_a_tmp*(1.0+redshift_241b)**ssfr_b_tmp)

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

# =============================================================================
# beta
# =============================================================================
#fig = plt.figure(figsize=(2*figuresize, figuresize))
fig = plt.figure(figsize=(4*figuresize, 2*figuresize))
ax1 = fig.add_axes([0, 1, 0.5, 0.5]) #[left, bottom, width, height]
ax1.plot(redshift, beta, label='This work')
ax1.fill_between(redshift, beta_16_arr, beta_84_arr, alpha=0.3)
#ax1.set_title('test')

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
ax1.plot(z_schreiber, beta_schreiber, label='Schreiber+15')
#ax1.plot(test_z, test_154_beta, label='154', marker='x', color='#1f77b4')
#ax1.plot(test_z, test_217_beta, label='217', marker='x', color='k')
#ax1.plot(test_z, test_154_beta_1, label='154', marker='v', color='#1f77b4', linestyle=':')
#ax1.plot(test_z, test_217_beta_1, label='217', marker='v', color='k', linestyle=':')
#ax1.plot(test_z, test_154_beta_2, label='154', marker='s', color='#1f77b4', linestyle=':')
#ax1.plot(test_z, test_217_beta_2, label='217', marker='s', color='k', linestyle=':')
#ax1.plot(redshift_217, beta_217, label='217', color='k')
#ax1.plot(redshift_tmp, beta_228, label='228')
#ax1.plot(redshift_tmp, beta_229, label='229 sfr cut, z1')
#ax1.plot(redshift_tmp, beta_230, label='230 z1')
#ax1.plot(redshift_tmp, beta_231, label='231 z1p5')
ax1.plot(redshift_tmp, beta_232, label='232')
ax1.plot(redshift_tmp, beta_233, label='233 sfr cut, z1')
ax1.plot(redshift_tmp, beta_234, label='234 z1')
ax1.plot(redshift_tmp, beta_235, label='235 z1p5')
ax1.plot(redshift_236, beta_236, label='236 0-1', color='r')
ax1.plot(redshift_237, beta_237, label='237 1-1p5', color='r')
ax1.plot(redshift_238, beta_238, label='238 1p5-2', color='r')
ax1.plot(redshift_239, beta_239, label='239 2-2p5', color='r')
ax1.plot(redshift_240, beta_240, label='240 2p5-3', color='r')
ax1.plot(redshift_241a, beta_241a, label='241a 3-6p5', color='r')
ax1.plot(redshift_241b, beta_241b, label='241b 3-6p5', color='r')

ax1.set_xlim(0.3, 7.7)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')

ax1.set_ylim(0.0, 2.0)
#ax1.set_ylim(0.3, 1.3)
ax1.set_ylabel(r'$\mathrm{MS\,Slope,}\,\beta$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

# =============================================================================
# alpha
# =============================================================================

ax2 = fig.add_axes([0, 0.5, 0.5, 0.5]) #[left, bottom, width, height]
ax2.plot(redshift, alpha, label='This work')
ax2.fill_between(redshift, alpha_16_arr, alpha_84_arr, alpha=0.3)


ax2.scatter(z_san, alpha_san_n, label='Santini+17')
ax2.errorbar(z_san, alpha_san_n, yerr=alpha_err_san_n, ls='none')
ax2.scatter(z_san0, alpha_san0_n, label='Santini+17 Raw')
ax2.errorbar(z_san0, alpha_san0_n, yerr=alpha_err_san0_n, ls='none')
ax2.scatter(z_salmon, alpha_salmon_n, label='Salmon+15')
ax2.scatter(z_steinhardt, alpha_steinhardt_n, label='Steinhardt+14')
ax2.scatter(z_kurc, alpha_kurc_n, label='Kurczynski+16')
ax2.plot(z_speagle, alpha_speagle_n, label='Speagle+14', linestyle=':')
ax2.plot(z_schreiber, alpha_schreiber_n, label='Schreiber+15')
#ax2.plot(test_z, test_154_9p7, label='154', marker='x', color='#1f77b4')
#ax2.plot(test_z, test_217_9p7, label='217', marker='+', color='k')
#ax2.plot(redshift_217, alpha_217, label='217', color='k')
#ax2.plot(redshift_tmp, alpha_228, label='228')
#ax2.plot(redshift_tmp, alpha_229, label='229 sfr cut, z1')
#ax2.plot(redshift_tmp, alpha_230, label='230 z1')
#ax2.plot(redshift_tmp, alpha_231, label='231 z1p5')
ax2.plot(redshift_tmp, alpha_232, label='232')
ax2.plot(redshift_tmp, alpha_233, label='233 sfr cut, z1')
ax2.plot(redshift_tmp, alpha_234, label='234 z1')
ax2.plot(redshift_tmp, alpha_235, label='235 z1p5')
ax2.plot(redshift_236, alpha_236, label='236 0-1', color='r')
ax2.plot(redshift_237, alpha_237, label='237 1-1p5', color='r')
ax2.plot(redshift_238, alpha_238, label='238 1p5-2', color='r')
ax2.plot(redshift_239, alpha_239, label='239 2-2p5', color='r')
ax2.plot(redshift_240, alpha_240, label='240 2p5-3', color='r')
ax2.plot(redshift_241a, alpha_241a, label='241a 3-6p5', color='r')
ax2.plot(redshift_241b, alpha_241b, label='241b 3-6p5', color='r')



p1 = np.random.uniform(size=1000)*10.0
p2 = np.random.uniform(size=1000)*10.0

for k in range(len(p1)):
    ax2.plot(redshift_tmp, np.log10(p1[k]*(1.0+redshift_tmp)**p2[k]) + 9.7 - 9.0, color='k', alpha=0.1)



ax2.set_xlim(0.3, 7.7)
ax2.set_xlabel('Redshift', labelpad=10)
ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax2.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax2.xaxis.set_tick_params(labelsize=fontsize_axes)

ax2.set_ylim(-0.7, 10)
ax2.set_ylabel(r'MS Normalisation, $\alpha$', labelpad=10)
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax2.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax2.yaxis.set_tick_params(labelsize=fontsize_axes)

ax2.legend(bbox_to_anchor=(1, 0), loc=4, frameon=False, fontsize=16)

# =============================================================================
# ssfr
# =============================================================================

ax3 = fig.add_axes([0, 0, 0.5, 0.5]) #[left, bottom, width, height]
ax3.plot(redshift, ssfr, label='This work')
ax3.fill_between(redshift, ssfr_16_arr, ssfr_84_arr, alpha=0.3)


ax3.plot(np.linspace(0.3, 7.7, 1000), np.log10((1+np.linspace(0, 7, 1000))**2.25))


ax3.set_xlim(0.3, 7.7)
ax3.set_xlabel('Redshift', labelpad=10)
ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax3.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax3.xaxis.set_tick_params(labelsize=fontsize_axes)

ax3.set_ylim(-2.5, 2.5)
ax3.set_ylabel(r'log(sSFR / Gyr)', labelpad=10)
ax3.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax3.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax3.yaxis.set_tick_params(labelsize=fontsize_axes)

ax3.legend(bbox_to_anchor=(1, 0), loc=4, frameon=False, fontsize=fontsize_legend)



# =============================================================================
# plt
# =============================================================================
plt.savefig('001_parameters_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

# Save figure
#plt.savefig('Final_Plot.png', dpi=300, transparent=False, bbox_inches='tight')


#%%
# =============================================================================
# redshift vs redshift inc IRAC flags
# =============================================================================
'''
#https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3
with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/recalc_chi2/astrodeep_pickle.p', 'rb') as f:
    astrodeep_pickle = pickle.load(f, encoding='latin1') 
#print(astrodeep_pickle.keys())
'''

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)

with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/data.p', 'rb') as f:
    data = pickle.load(f, encoding='latin1') 
#print(data.keys())

# =============================================================================
# subsets
# =============================================================================

idx1 = (data['field_AD']%2.0==0.0) # clusters
idx2 = (data['relflag_AD']==1.0) # relflag
idx3 = (-2.5*np.log10(AD['b_H160']*1e-6) + 8.90 < 27.5) # H band cut

idx3_IRAC = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)
idx3_IRAC = np.logical_and(idx3_IRAC, data['redshift_BEAGLE']>4.0)
idx3_IRAC = ~idx3_IRAC

idx3_IRAC_removed = np.logical_and(data['CH1_BEAGLE_input']<-60.0, data['CH2_BEAGLE_input']<-60.0)

MCLmassLow = np.empty(data['redshift_BEAGLE'].shape)
for i in range(len(MCLmassLow)):
    if data['redshift_BEAGLE'][i] <2.1789654:
        MCLmassLow[i] = 8.0
    elif data['redshift_BEAGLE'][i] > 4.195:
        MCLmassLow[i] = 9.0
    else:
        MCLmassLow[i] = 6.91926521 + 0.49598529*data['redshift_BEAGLE'][i]
idx4 = (data['mass_BEAGLE_stellar'] + data['mag_AD'] > MCLmassLow) 

idx5_1 = (abs(data['redshift_BEAGLE']-((0.5+6.5)/2.0)) < (((0.5+6.5)/2.0) - 0.5)) # 0.5<redshift_BEAGLE<8.0 (only affects up to visual inspection of z=3.5)
idx5_2 = (abs(data['redshift_BEAGLE']-data['redshift_AD']) < 1.0) # |redshift_BEAGLE-redshift_AD|<1
idx5_3 = np.logical_and((data['redshift_BEAGLE'] < 3.5), (data['redshift_AD'] < 3.5)) 
idx5_zlt3p5 = np.logical_and(idx5_1,idx5_2)
idx5_zlt3p5 = np.logical_and(idx5_zlt3p5,idx5_3)
vis = np.genfromtxt('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/redshift_investigation_4/investigation_3_and_4_selection.csv', delimiter=",", skip_header=1)   
idx5_zgt3p5 = np.full(len(data['id_AD']), False)
for i in range(len(data['id_AD'])):
    idx5_zgt3p5_temp = np.isin(data['id_AD'][i],vis[:,1][(vis[:,0]==data['field_AD'][i])&(vis[:,3]==1)])
    if idx5_zgt3p5_temp:
        idx5_zgt3p5[i] = True   
idx5_zgt3p5 = np.logical_and(idx5_zgt3p5, data['redshift_BEAGLE'] < 6.5)
idx5_z = np.logical_or(idx5_zlt3p5,idx5_zgt3p5)
   
TmassHigh = 9.244 + (0.753*np.minimum(data['redshift_BEAGLE'], 4.0)) - (0.090*(np.minimum(data['redshift_BEAGLE'], 4.0)**2)) # all galaxies (used by santini)
idx6 = (data['mass_BEAGLE_stellar'] < (TmassHigh))

idx7 = (data['new_min_chi2_BEAGLE']>0) & (data['new_min_chi2_BEAGLE']<13.28) # chi-squared
idx8 = (data['sfr_BEAGLE_instant']>-5.0) & (data['sfr_BEAGLE_instant']<10.0) # arbitrary upper and lower sfr cuts
idx9 = (data['xsig_GMM_3d'][:,0]>-100.0) # GMM 3d not fitted are assigned x y z sig of -103

# =============================================================================
# plot - point of this is to compare photo-z, so need it to be on objects we would expect a decent photo z
'''
can include parallels and clusters
probably do H cut, chi2 cut, THIS IS WITH ME to work out what should be included

anything here needs to be independent of visual inspection except for final yellows

be nice to get the purple ones as a contour plot, check for the one I can imagine

'''
# =============================================================================


# relflag + H band cut
# chi2

# ignore parallel cut
# ignore IRAC cut
# ignore mass cut
# ignore redshift cuts
# ignore arbitrary sfr cut
# ignore GMM cut

idx = np.logical_and(idx2,idx3) 
idx = np.logical_and(idx,idx7)


###########
nbins_x = 30
nbins_y = 30
min_x = np.min(data['redshift_AD'][idx])
max_x = np.max(data['redshift_AD'][idx])
binsize_x = (max_x-min_x)/nbins_x
min_y = np.min(data['redshift_BEAGLE'][idx])
max_y = np.max(data['redshift_BEAGLE'][idx])
binsize_y = (min_y-max_y)/nbins_y
data1 = plt.hist2d(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], bins=[40,40], range=[[min_x,max_x],[min_y,max_y]])

###########


fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
      
ax1.contourf((data1[1][1:]+data1[1][:-1])/2., \
                (data1[2][1:]+data1[2][:-1])/2., np.log10(np.transpose(data1[0])), \
                 cmap=cm.gist_yarg)
#ax1.hist2d(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], bins=60, range=[[0,9.8],[0,9.8]], norm=mpl.colors.LogNorm(), cmap=mpl.cm.binary)
#ax1.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], alpha=0.3)

# ~IRAC cut

idx_IRAC = np.logical_and(idx,idx3_IRAC_removed)
ax1.scatter(data['redshift_AD'][idx_IRAC], data['redshift_BEAGLE'][idx_IRAC], s=1, alpha=1.0, color='blue')

# IRAC cut
# parallel cut
# mass cut
# redshift cuts + visual inspection
# arbitrary sfr cut
# GMM cut

idx = np.logical_and(idx,idx3_IRAC)
idx = np.logical_and(idx,idx1)
idx = np.logical_and(idx,idx4)
idx = np.logical_and(idx,idx6)
idx = np.logical_and(idx,idx5_z)
idx = np.logical_and(idx,idx8)
idx = np.logical_and(idx,idx9)

ax1.scatter(data['redshift_AD'][idx], data['redshift_BEAGLE'][idx], s=1, alpha=1.0, color='red')

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

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)
plt.savefig('002_redshift_vs_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()








# =============================================================================
# Boogard #1
# =============================================================================
#%%
with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
    data = pickle.load(f, encoding='latin1') 
print(data.keys())

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((9.7, 9.7), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

x_test = np.linspace(5, 12, 10)
#z = 1
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)
#z = 2
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)
#z = 3
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)
#z = 4
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)
#z = 5
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)
#z = 6
#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='r', alpha=0.5)


idx_sort = np.argsort(data['redshift_BEAGLE'])
scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], data['sfr_BEAGLE_instant'][idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=0, vmax=7, alpha=1.0)

#ax1.set_title('test')

ax1.set_xlim(6.5, 10.5)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)


cb.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label('Redshift', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)
plt.savefig('003_MS_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

# =============================================================================
# Boogard #2
# =============================================================================
#%%

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]


#ax1.set_title('test')

sfr_surface = ((beta_a + beta_b*data['redshift_BEAGLE'])*(data['mass_BEAGLE_stellar']-9.7)) + np.log10(ssfr_a*(1.0+data['redshift_BEAGLE'])**ssfr_b) + 9.7 - 9.0 
idx_sort = np.argsort(data['redshift_BEAGLE'])
#scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-sfr_surface)[idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=0, vmax=7, alpha=0.6)


z = 1

sfr_surface_norm = ((beta_a+z*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b))
sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+9.7-9.0-(9.7*(beta_a+data['redshift_BEAGLE']*beta_b))
offset = sfr_surface_real - sfr_surface_norm
#scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-offset)[idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=0, vmax=7, alpha=1.0)
scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=0, vmax=7, alpha=1.0)


ax1.plot(x_test, x_test*0, color='k', alpha=0.5)

ax1.set_xlim(5.5, 11.5)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(r'$\Delta_{MS}$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)


cb.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label('Redshift', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

plt.savefig('004_MS_collapsed_redshift.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


# =============================================================================
# Boogard #5 colourcode for Hogg
# =============================================================================
#%%

sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+9.7-9.0-(9.7*(beta_a+data['redshift_BEAGLE']*beta_b))

pbad          = 0.15
outlier_mean  = 1.0
outlier_sigma = 1.0
sig0          = 0.3

log_p_xi_eta_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
p_bad = norm.pdf(data['sfr_BEAGLE_instant'][idx_sort], scale=outlier_sigma, loc=outlier_mean)

z_bad = pbad*p_bad
z_good = (1-pbad)*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
    data = pickle.load(f, encoding='latin1') 
print(data.keys())

fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

ax1.plot((0, 20), (-10, 10), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((0, 20), (-9, 11), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((0, 20), (-8, 12), color='k', alpha=0.5, linestyle='dashed')
ax1.plot((9.7, 9.7), (-10, 10), color='k', alpha=0.5, linestyle='dotted')

scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], data['sfr_BEAGLE_instant'][idx_sort], c=np.log10(z_good[idx_sort]/z_bad[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-0.5, vmax=2.5)

ax1.set_xlim(6.5, 10.5)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(r'$\log(\psi \, / \, \mathrm{M_{\odot} \, yr^{-1}})$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)


cb.set_ticks(np.linspace(0, 2, 3))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

plt.savefig('005_MS_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()

#plt.scatter(z_good, z_bad)
#plt.hist(p_bad)
#plt.hist(np.exp(log_p_eta_xi_theta), bins=100)
#plt.hist(log_p_eta_xi_theta, bins=100)
#print(max(log_p_eta_xi_theta))
#print(z_good/z_bad)
#plt.hist(z_good/z_bad, bins=10000)
#plt.xlim(0, 2)
#plt.show()
#
#plt.hist(z_good-z_bad, bins=10000)
#plt.xlim(0, 2)
#plt.show()

# =============================================================================
# Boogard #6 colourcode for Hogg
# =============================================================================
#%%


sfr_surface_real = ((beta_a+data['redshift_BEAGLE']*beta_b)*data['mass_BEAGLE_stellar'])+(np.log10(ssfr_a*((1+data['redshift_BEAGLE'])**ssfr_b)))+9.7-9.0-(9.7*(beta_a+data['redshift_BEAGLE']*beta_b))

pbad          = 0.15
outlier_mean  = 1.0
outlier_sigma = 1.0
sig0          = 0.3

log_p_xi_eta_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
log_p_eta_xi_theta = norm.logpdf(data['sfr_BEAGLE_instant'], scale=sig0, loc=sfr_surface_real)
p_bad = norm.pdf(data['sfr_BEAGLE_instant'][idx_sort], scale=outlier_sigma, loc=outlier_mean)

z_bad = pbad*p_bad
z_good = (1-pbad)*np.exp(log_p_eta_xi_theta)

idx_sort = np.argsort(z_good/z_bad)

with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.p', 'rb') as f:
    data = pickle.load(f, encoding='latin1') 
#print(data.keys())




fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]

z = 1
beta_adj = (beta_a+z*beta_b) - (beta_a+data['redshift_BEAGLE']*beta_b)
sfr_adj = np.log10((10**data['sfr_BEAGLE_instant'])*((10**data['mass_BEAGLE_stellar'])**beta_adj)*(((1+z)/(1+data['redshift_BEAGLE']))**ssfr_b)) - 9.7*beta_adj


#scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], sfr_adj[idx_sort], c=data['redshift_BEAGLE'][idx_sort], vmin=0, vmax=7, alpha=1.0)

#scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], sfr_adj[idx_sort], c=z_good[idx_sort]/z_bad[idx_sort], cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=0, vmax=2)

scatter = ax1.scatter(data['mass_BEAGLE_stellar'][idx_sort], (data['sfr_BEAGLE_instant']-sfr_surface_real)[idx_sort], c=np.log10(z_good[idx_sort]/z_bad[idx_sort]), cmap=mpl.cm.get_cmap('coolwarm'), alpha=1.0, vmin=-0.5, vmax=2.5)


#plt.hist(np.log10(z_good[idx_sort]/z_bad[idx_sort]), bins=100, range=[0, 10])
#plt.show()


#ax1.plot(x_test, ((beta_a+z*beta_b)*x_test)+(np.log10(ssfr_a*((1+z)**ssfr_b)))+9.7-9.0-(9.7*(beta_a+z*beta_b)), color='k', alpha=0.5)
ax1.plot(x_test, x_test*0, color='k', alpha=0.5)

ax1.set_xlim(5.5, 11.5)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
ax1.set_ylabel(r'$\Delta_{MS}$', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

#ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

cbaxes = fig.add_axes([0.85, 0, 0.05, 0.84]) #[left, bottom, width, height]
cb = plt.colorbar(scatter, cax = cbaxes)


cb.set_ticks(np.linspace(0, 2, 3))
#cb.set_ticklabels(['', 1, 2, 3, 4, 5, 6, ''])
cb.set_label(r'$\log(\mathrm{p}_{MS} \, / \, \mathrm{p}_{OL})$', rotation=270, labelpad=30)
cb.ax.tick_params(axis='y', size=5, width=2, direction='in', labelsize=fontsize_axes)

plt.savefig('006_MS_collapsed_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()




# =============================================================================
# bias test histograms
# =============================================================================
#%%

tests = ['110', '111', '112', '113', '114', '115', '116', '117', '118', '119']
tests = ['115', '117', '110', '119', '112', '111', '116', '113', '118', '114']
names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
names_plot = ['alphaN a', 'alphaN b', 'beta a', 'beta b', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']
trues = [4.0, -0.75, 0.5, 0.25, 0.3, 0.3, 0.0, 2.0] # outlier_mean is REALLY 0
'''
dic = {}
for name in names:
    dic[name+'_16'] = []
    dic[name] = []
    dic[name+'_84'] = []
    
for test in tests:

    # python 3
    with open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_truncated_{}.p'.format(test), 'rb') as f:
        chain_bias = pickle.load(f, encoding='latin1') 
    
    nChains=4
    minIter=len(chain_bias)/nChains
    burn=int(0.75*minIter)
#        burn=0
    chain_arr = []
    for i in range(nChains):
        start = int(minIter*i+burn)
        finish = int(minIter*(i+1))
#        print(start, finish)
        chain_arr.append(chain_bias[start:finish])
    chain_bias = np.concatenate(chain_arr)   

    for name in names:
        dic[name+'_16'].append(float('{:.3f}'.format(np.percentile(chain_bias[name], 16))))
        dic[name].append(float('{:.3f}'.format(np.median(chain_bias[name]))))
        dic[name+'_84'].append(float('{:.3f}'.format(np.percentile(chain_bias[name], 84))))
'''
#%%
#plot absolute around 0
x_values = np.linspace(1.5, 2., 10)
fig = plt.figure(figsize=(2*figuresize, 0.5*figuresize))
ax1 = fig.add_axes([0, 0, 1, 1]) #[left, bottom, width, height]
for j, name in enumerate(names):
    ax1.scatter(x_values+j, np.array(dic[name]) - trues[j], color='k', marker='x')
    ax1.plot((x_values+j, x_values+j), (np.array(dic[name+'_16']) - trues[j], np.array(dic[name+'_84']) - trues[j]), color='k', zorder=0, alpha=0.5)   
ax1.plot((0,10),(0,0), color='k', alpha=0.5)

ax1.set_xlim(1.25, 9.25)
#ax1.set_xlabel(r'TEST', labelpad=10)
#plt.xticks(np.linspace(1.75, 8.75, len(names_plot)), names_plot)
ax1.set_xticks(np.linspace(1.75, 8.75, len(names_plot)))
ax1.set_xticklabels(names_plot)
#ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
#ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
#ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-0.75, 0.75)
ax1.set_ylabel(r'TEST', labelpad=10)
ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
ax1.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', left='on', right='on')
ax1.yaxis.set_tick_params(labelsize=fontsize_axes)

plt.savefig('007_bias_test_histograms.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#
## =============================================================================
## Corner Plots
## =============================================================================
#%%

names = ['alphaN_a', 'alphaN_b', 'beta_a', 'beta_b', 'sig0', 'pbad', 'outlier_mean', 'outlier_sigma']
names_plot = ['alphaN a', 'alphaN b', 'beta a', 'beta b', 'sig0', 'pbad', 'outlier mean', 'outlier sigma']

'''
# python 3
#with open('/Users/lester/Documents/linmix_files/lm_chain_mock_hogg_redshift_truncated_{}.p'.format('110'), 'rb') as f:
#    chain_corner = pickle.load(f, encoding='latin1') 

#nChains=4
#minIter=len(chain_corner)/nChains
#burn=int(0.75*minIter)
##        burn=0
#chain_arr = []
#for i in range(nChains):
#    start = int(minIter*i+burn)
#    finish = int(minIter*(i+1))
##        print(start, finish)
#    chain_arr.append(chain_corner[start:finish])
#chain_corner = np.concatenate(chain_arr)   

with open('/Users/lester/Documents/linmix_files/lm_chain_3d_Hogg_truncated_normal_sSFR_delayed_uniform_logtau_scenario_23_data_z0p5_154_subset_002.p', 'rb') as f:
    chain_corner = pickle.load(f, encoding='latin1') 
'''

data = np.array([chain_corner['alphaN_a'],chain_corner['alphaN_b'],chain_corner['beta_a'],chain_corner['beta_b'],chain_corner['sig0'],chain_corner['pbad'],chain_corner['outlier_mean'],chain_corner['outlier_sigma']]).T

# Plot it.
#fig1 = plt.figure(figsize=(figuresize, figuresize))
#fig1 = plt.figure()
#figure = corner.corner(data, labels=names_plot,
#                       quantiles=[0.16, 0.5, 0.84],
#                       show_titles=True, title_kwargs={"fontsize": fontsize_axes}, fig=fig1)

figure = corner.corner(data, labels=names_plot,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": fontsize_axes})
plt.savefig('008_corner_plots.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#        chain_burn = pickle.load(open(,'r'))


# =============================================================================
# heatplots from GMM
# =============================================================================
#%% 

#154
#ssfr_a_hp = 0.056
#ssfr_b_hp = 2.54
#beta_a_hp = 0.687
#beta_b_hp = -0.014

##218 1st chain 4000burn
#ssfr_a_hp = 0.003
#ssfr_b_hp = 6.831
#beta_a_hp = 0.446
#beta_b_hp = 0.234
#
##219 1st chain 4000burn median mass > 8.0
#ssfr_a_hp = 0.004
#ssfr_b_hp = 6.44
#beta_a_hp = 0.629
#beta_b_hp = 0.139
#
##220 1st chain 4000burn median mass > 8.5
#ssfr_a_hp = 0.006
#ssfr_b_hp = 6.435
#beta_a_hp = 0.295
#beta_b_hp = 0.589
#
##221 1st chain 4000burn median mass > 9.0
#ssfr_a_hp = 0.001
#ssfr_b_hp = 8.583
#beta_a_hp = -1.514
#beta_b_hp = 2.011
#
#217 1st chain 4000burn -2<sfr<2
#ssfr_a_hp = 0.07
#ssfr_b_hp = 2.351
#beta_a_hp = 0.712
#beta_b_hp = -0.075


num = 3
z_med_hp = 1.25
z_med_hp_gap = 0.25
santini_idx = int(z_med_hp-z_med_hp_gap) # 1,2,3,4,5

#alpha_hp = np.log10(ssfr_a_hp*(1.0+z_med_hp)**ssfr_b_hp) + 9.7 - 9.0
#beta_hp = beta_a_hp + z_med_hp*beta_b_hp
#z_med_hpp = z_med_hp + z_med_hp_gap
#alpha_hpp = np.log10(ssfr_a_hp*(1.0+z_med_hpp)**ssfr_b_hp) + 9.7 - 9.0
#beta_hpp = beta_a_hp + z_med_hpp*beta_b_hp
#z_med_hpm = z_med_hp - z_med_hp_gap
#alpha_hpm = np.log10(ssfr_a_hp*(1.0+z_med_hpm)**ssfr_b_hp) + 9.7 - 9.0
#beta_hpm = beta_a_hp + z_med_hpm*beta_b_hp


#154 original scenario 23 with Hogg etc
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.056,2.54,0.687,-0.014
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_154 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_154 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_154 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#217 2x 3750 to 5000 without Hogg across all redshifts, -2<sfr<2
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.07,2.351,0.712,-0.075
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_217 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_217 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_217 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)
 
#232 2x 3750 to 5000 ALL
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.048,2.682,0.659,-0.003
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_232 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_232 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_232 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#233 2x 3750 to 5000 sfr cut, z>1
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.137,1.854,0.766,-0.019
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_233 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_233 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_233 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#234 2x 3750 to 5000 z>1
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.164,1.692,0.837,-0.053
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_234 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_234 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_234 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)

#235 2x 3750 to 5000 z>1.5
ssfr_a_tmp,ssfr_b_tmp,beta_a_tmp,beta_b_tmp = 0.315,1.209,0.975,-0.091
redshift_tmp = np.linspace(z_med_hp-z_med_hp_gap, z_med_hp+z_med_hp_gap, num)
beta_235 = beta_a_tmp + redshift_tmp*beta_b_tmp
alpha_235 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp) + 9.7 - 9.0
ssfr_235 = np.log10(ssfr_a_tmp*(1.0+redshift_tmp)**ssfr_b_tmp)


fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.fits'
#s23 = fits.open(fileName)
#print(s23.info())
#print(s23[1].header)
s23 = fits.open(fileName)[1].data # scenario 23

#idx_rdm = np.random.randint(0, len(s23)-1, 100) # random 20 objects
#idx_rdm = np.arange(len(s23))[abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5)&(s23['mass_BEAGLE_stellar']>8.0)] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5)&(s23['mass_BEAGLE_stellar']>8.5)] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5)&(s23['mass_BEAGLE_stellar']>9.0)] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['sfr_BEAGLE_instant'])<2.0)] # sfr cut
idx_rdm = np.arange(len(s23))[:] # all objects

x_hp = np.array([])
y_hp = np.array([])
z_hp = np.array([])

n_hp = 300 # number of samples to take from GMM in total
#for i in range(len(s23['id_GMM_3d'])):
for i in idx_rdm:
    G_arr = np.random.choice([0,1,2], 10, p=s23['amp_GMM_3d'][i])

    for G in range(3):
        
        mean = np.array([s23['x_GMM_3d'][i,G],s23['y_GMM_3d'][i,G],s23['z_GMM_3d'][i,G]])
        cov = np.array([[np.power(s23['xsig_GMM_3d'][i,G],2), s23['xycov_GMM_3d'][i,G], s23['xzcov_GMM_3d'][i,G]],[s23['xycov_GMM_3d'][i,G], np.power(s23['ysig_GMM_3d'][i,G],2), s23['yzcov_GMM_3d'][i,G]],[s23['xzcov_GMM_3d'][i,G], s23['yzcov_GMM_3d'][i,G], np.power(s23['zsig_GMM_3d'][i,G],2)]])

        xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s23['amp_GMM_3d'][i,G]))

        x_hp = np.concatenate((x_hp,xyz[:,0]))
        y_hp = np.concatenate((y_hp,xyz[:,1]))
        z_hp = np.concatenate((z_hp,xyz[:,2]))

#plt.figure(figsize=(10, 6))
#xlow = 5
#xhigh = 12
#ylow = -5
#yhigh = 5
#plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
#plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
#plt.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
#plt.colorbar()
#cmap = cm.get_cmap('viridis')
#rgba = cmap(0)
#ax = plt.axes()
#ax.set_facecolor(rgba)
#plt.scatter(s23['mass_BEAGLE_stellar'][idx_rdm], s23['sfr_BEAGLE_instant'][idx_rdm], color='r', marker='x')
#plt.show()



x_hp = x_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
y_hp = y_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]
z_hp = z_hp[abs(z_hp - z_med_hp) < z_med_hp_gap]


### litte bit to find out the value I think I should pass through at mass 9.7
#
#y_7p7 = y_hp[abs(x_hp-7.7)<0.1]
#y_8p7 = y_hp[abs(x_hp-8.7)<0.1]
#y_9p7 = y_hp[abs(x_hp-9.7)<0.1]
##print(np.mean(y_7p7),np.mean(y_9p7))
#print(np.mean(y_8p7))
### ### ###


fig = plt.figure(figsize=(2*figuresize, 2*figuresize))

ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
ax1.set_title('{} - {}'.format(z_med_hp-
z_med_hp_gap, z_med_hp+
z_med_hp_gap))

xlow = 5.5
xhigh = 11.5
ylow = -8.5
yhigh = 3.5
h = ax1.hist2d(x_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
ax1.scatter(s23['mass_BEAGLE_stellar'][idx_rdm][abs(s23['redshift_BEAGLE'] - z_med_hp) < z_med_hp_gap], s23['sfr_BEAGLE_instant'][idx_rdm][abs(s23['redshift_BEAGLE'] - z_med_hp) < z_med_hp_gap], color='r', marker='x')
ax1.plot((xlow,xhigh), (alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xlow,alpha_san[santini_idx-1] + beta_san[santini_idx-1]*xhigh), color='w') # santini
#ax1.plot((-0.3,9.7), (-6.5,0.5)) # my work 1<z<2
#ax1.plot((-0.3,9.7), (-9.5,1.0)) # santini 1<z<2
#ax1.plot((-0.3,9.7), (-4.7,1.3)) # my work 4<z<5
#ax1.plot((-0.3,9.7), (-8.2,1.3)) # santini 4<z<5
#ax1.plot((7.7,8.7,9.7), (test_154_7p7[int(z_med_hp-z_med_hp_gap)],test_154_8p7[int(z_med_hp-z_med_hp_gap)],test_154_9p7[int(z_med_hp-z_med_hp_gap)]), color='w', marker='v')
#ax1.plot((7.7,8.7,9.7), (test_217_7p7[int(z_med_hp-z_med_hp_gap)],test_217_8p7[int(z_med_hp-z_med_hp_gap)],test_217_9p7[int(z_med_hp-z_med_hp_gap)]), color='w', marker='s') 
ax1.plot((9.7,9.7), (ylow, yhigh), color='w') 
#ax1.plot((xlow,xhigh), (alpha_hp + beta_hp*(xlow-9.7),alpha_hp + beta_hp*(xhigh-9.7)), color='r') # my work
#ax1.plot((xlow,xhigh), (alpha_hpp + beta_hpp*(xlow-9.7),alpha_hpp + beta_hpp*(xhigh-9.7)), color='r', linestyle=':') # my work
#ax1.plot((xlow,xhigh), (alpha_hpm + beta_hpm*(xlow-9.7),alpha_hpm + beta_hpm*(xhigh-9.7)), color='r', linestyle=':') # my work

ax1.plot((xlow,xhigh), (alpha_154[1] + beta_154[1]*(xlow-9.7),alpha_154[1] + beta_154[1]*(xhigh-9.7)), color='r', label='154 original') # my work
ax1.plot((xlow,xhigh), (alpha_217[1] + beta_217[1]*(xlow-9.7),alpha_217[1] + beta_217[1]*(xhigh-9.7)), color='k', label='217 no Hogg, sfr cut')
ax1.plot((xlow,xhigh), (alpha_232[1] + beta_232[1]*(xlow-9.7),alpha_232[1] + beta_232[1]*(xhigh-9.7)), linestyle='dashed', label='232')
ax1.plot((xlow,xhigh), (alpha_233[1] + beta_233[1]*(xlow-9.7),alpha_233[1] + beta_233[1]*(xhigh-9.7)), linestyle='dashed', label='233 z1, sfr cut')
ax1.plot((xlow,xhigh), (alpha_234[1] + beta_234[1]*(xlow-9.7),alpha_234[1] + beta_234[1]*(xhigh-9.7)), linestyle='dashed', label='234 z1')
ax1.plot((xlow,xhigh), (alpha_235[1] + beta_235[1]*(xlow-9.7),alpha_235[1] + beta_235[1]*(xhigh-9.7)), linestyle='dashed', label='235 z1p5')

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

ax1.legend(bbox_to_anchor=(1, 1), loc=1, frameon=False, fontsize=16)

plt.savefig('00006_MS_collapsed_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#%%
fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
xlow = 5.5
xhigh = 11.5
ylow = 0
yhigh = 7
h = ax1.hist2d(x_hp, z_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
ax1.scatter(s23['mass_BEAGLE_stellar'][idx_rdm], s23['redshift_BEAGLE'][idx_rdm], color='r', marker='x')
ax1.set_xlim(5.5, 11.5)
ax1.set_xlabel(r'$\log(m_{\star} \, / \, \mathrm{M_{\odot}})$', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(0, 7)
ax1.set_ylabel(r'Redshift', labelpad=10)
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

plt.savefig('00006_MS_collapsed_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()



fig = plt.figure(figsize=(figuresize, figuresize))
ax1 = fig.add_axes([0, 0, 0.85, 0.84]) #[left, bottom, width, height]
xlow = 0
xhigh = 7
ylow = -3.5
yhigh = 3.5
h = ax1.hist2d(z_hp, y_hp, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
ax1.scatter(s23['redshift_BEAGLE'][idx_rdm], s23['sfr_BEAGLE_instant'][idx_rdm], color='r', marker='x')
ax1.set_xlim(0, 7)
ax1.set_xlabel(r'Redshift', labelpad=10)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax1.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', bottom='on', top='on')
ax1.xaxis.set_tick_params(labelsize=fontsize_axes)

ax1.set_ylim(-3.5, 3.5)
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

plt.savefig('00006_MS_collapsed_hogg.png', dpi=300, transparent=False, bbox_inches='tight')
plt.show()


#%%
# =============================================================================
# GMM individuals
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/kelly/delayed_uniform_logtau/scenario_23_data_z0p5.fits'
#s23 = fits.open(fileName)
#print(s23.info())
#print(s23[1].header)
s23 = fits.open(fileName)[1].data # scenario 23

#idx_rdm = np.random.randint(0, len(s23)-1, 100) # random 20 objects
#idx_rdm = np.arange(len(s23))[abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5)&(s23['mass_BEAGLE_stellar']>8.0)] # redshift bin of objects
idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - 0.5) < 0.5)&(s23['mass_BEAGLE_stellar']>8.5)] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['redshift_BEAGLE'] - z_med_hp) < 0.5)&(s23['mass_BEAGLE_stellar']>9.0)] # redshift bin of objects
#idx_rdm = np.arange(len(s23))[(abs(s23['sfr_BEAGLE_instant'])<2.0)] # sfr cut
#idx_rdm = np.arange(len(s23))[:] # all objects



n_hp = 500 # number of samples to take from GMM in total
#for i in range(len(s23['id_GMM_3d'])):
for i in idx_rdm:
#for i in range(4):
    G_arr = np.random.choice([0,1,2], 10, p=s23['amp_GMM_3d'][i])

    x_hp = np.array([])
    y_hp = np.array([])
    z_hp = np.array([])

    for G in range(3):
        
        mean = np.array([s23['x_GMM_3d'][i,G],s23['y_GMM_3d'][i,G],s23['z_GMM_3d'][i,G]])
        cov = np.array([[np.power(s23['xsig_GMM_3d'][i,G],2), s23['xycov_GMM_3d'][i,G], s23['xzcov_GMM_3d'][i,G]],[s23['xycov_GMM_3d'][i,G], np.power(s23['ysig_GMM_3d'][i,G],2), s23['yzcov_GMM_3d'][i,G]],[s23['xzcov_GMM_3d'][i,G], s23['yzcov_GMM_3d'][i,G], np.power(s23['zsig_GMM_3d'][i,G],2)]])

        xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s23['amp_GMM_3d'][i,G]))

        x_hp = np.concatenate((x_hp,xyz[:,0]))
        y_hp = np.concatenate((y_hp,xyz[:,1]))
        z_hp = np.concatenate((z_hp,xyz[:,2]))
        
    plt.title(int(s23['id_GMM_3d'][i]))
    plt.hist(x_hp-9.0, histtype=u'step', density=True, label='mass - 9')
    plt.hist(y_hp, histtype=u'step', density=True, label='sfr')
    plt.hist(z_hp, histtype=u'step', density=True, label='z')
#    plt.legend()
    plt.show()




