#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# =============================================================================
# NOTES
# =============================================================================

'''
NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

THESE BECOME NAN: BEAGLE had -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted

# log(0) -> -inf (mass_AD, mass_AD_neb, sfr_AD, sfr_AD_neb)
# lof(-ve) -> nan (mass_BEAGLE_tot and )

'''

# =============================================================================
# OPTIONS
# =============================================================================

# Clusters, H<27.5, 1.3<redshift_AD<2.0, mass_AD_neb>8.4, relflag==1, mass_AD_neb<10.2

# SAME + 1.3<redshift_BEAGLE<2.0
option_BEAGLE_redshift = 0

# SAME + mass_BEAGLE_stellar
option_BEAGLE_mass = 0

# Changed limit on BEAGLE mass to 7.863
option_BEAGLE_mass_lower = 0

# Adjusting upper limit by same factor (8.4 - 7.863)
option_BEAGLE_mass_upper = 0

# SAME + sfr_BEAGLE_instant
option_BEAGLE_sfr = 1

# all fields
option_all_fields = 1 # otherwise default is clusters

# =============================================================================
# # Santini approx 1711 sources 1.3 < z < 6.0
# =============================================================================

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)

sbf = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/npy_files_matching_AD_and_BEAGLE_8_fields/'
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'

with np.errstate(divide='ignore', invalid='ignore'):
    data = {    'field_AD':             AD['field'],
                'id_AD':                AD['ID'], 
                'mag_AD':               np.log10(AD['MAGNIF']), 
                'redshift_AD':          AD['ZBEST'], 
                'mass_AD':              np.log10(AD['MSTAR']*1e9) - np.log10(AD['MAGNIF']), 
                'mass_AD_neb':          np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']), 
                'sfr_AD':               np.log10(AD['SFR']) - np.log10(AD['MAGNIF']), 
                'sfr_AD_neb':           np.log10(AD['SFR_NEB']) - np.log10(AD['MAGNIF']), 
                'relflag_AD':           AD['RELFLAG'],
                'id_BEAGLE':            np.load(sbf+'id_BEAGLE.npy'), 
                'sfr_SAN':              np.load(sfr_SAN_location) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_tot':      np.log10(np.load(sbf+'mass_BEAGLE_tot.npy')) - np.log10(AD['MAGNIF']), 
                'mass_BEAGLE_stellar':  np.log10(np.load(sbf+'mass_BEAGLE_stellar.npy')) - np.log10(AD['MAGNIF']), 
                'sfr_BEAGLE_instant':   np.load(sbf+'sfr_BEAGLE_instant.npy') - np.log10(AD['MAGNIF']), 
                'redshift_BEAGLE':      np.load(sbf+'redshift_BEAGLE.npy'), 
                'tau_BEAGLE':           np.load(sbf+'tau_BEAGLE.npy'), 
                'tauv_BEAGLE':          np.load(sbf+'tauv_BEAGLE.npy'), 
                'msa_BEAGLE':           np.load(sbf+'msa_BEAGLE.npy'), 
                'metallicity_BEAGLE':   np.load(sbf+'metallicity_BEAGLE.npy')                
    
                }

# =============================================================================
# clusters only?
# =============================================================================
if option_all_fields == 1:
    idx1 = (AD['field']%1.0==0.0) # all
elif option_all_fields == 0:
    idx1 = (AD['field']%2.0==0.0) # clusters


# =============================================================================
# H band mag cut of < 27.5
# =============================================================================
idx2 = (AD['H160']<27.5)
idx = np.logical_and(idx1,idx2)

# =============================================================================
# Photo-z from median of 6 methods
# =============================================================================
zLow = 1.3
zHigh = 2.0
wLim = (zHigh - zLow) / 2.0
zLim = zLow + wLim
idx3 = (abs(AD['ZBEST']-zLim) < wLim)
idx = np.logical_and(idx,idx3)

# =============================================================================
# Need sample to be complete above given mass (magnification NOT already included)
# =============================================================================
if option_BEAGLE_mass == 0:
    idx4 = (data['mass_AD_neb'] + np.log10(AD['MAGNIF']) > 8.4)
    idx = np.logical_and(idx,idx4)
elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_lower== 0:
    # changing mass to BEAGLE
    idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > 8.4)
    idx = np.logical_and(idx,idx4)
elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_lower == 1:
    # adjusting for offset in mass between AD and BEAGLE to get same number of objects (404)
    idx4 = (data['mass_BEAGLE_stellar'] + np.log10(AD['MAGNIF']) > 7.863)
    idx = np.logical_and(idx,idx4)

# =============================================================================
# Adding RELFLAG just to be sure...
# =============================================================================
idx5 = (AD['RELFLAG']==1.0)
idx = np.logical_and(idx,idx5)

# =============================================================================
# Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
# =============================================================================
if option_BEAGLE_mass == 0:
    idx6 = (np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']) < 10.2)
    idx = np.logical_and(idx,idx6)
elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_upper == 0:
    # changing mass to BEAGLE
    idx6 = (data['mass_BEAGLE_stellar'] < 10.2)
    idx = np.logical_and(idx,idx6)    
elif option_BEAGLE_mass == 1 and option_BEAGLE_mass_upper == 1:
    # adjusting for offset in mass between AD and BEAGLE 
    idx6 = (data['mass_BEAGLE_stellar'] < 10.2-(8.4-7.863))
    idx = np.logical_and(idx,idx6)

# =============================================================================
# BEAGLE REDSHIFT AS WELL AS ASTRODEEP
# =============================================================================
if option_BEAGLE_redshift == 1:
    idx7 = (abs(data['redshift_BEAGLE']-zLim) < wLim)
    idx = np.logical_and(idx,idx7)

# =============================================================================
# COMBINE THE ABOVE
# =============================================================================

#for key in data.keys():
#    print(key, len(data[key]), min(data[key]))
#    print(data[key])

for key in data.keys():
    data[key] = data[key][idx]
#    print(key, len(data[key]), min(data[key]))

# =============================================================================
# Santini
# =============================================================================
# logSFR = alpha log(M / M_9.7) + beta
z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

alpha_san = B_san - 9.7*A_san
beta_san = A_san
# CHECK THIS
alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_err_san = A_err_san


# =============================================================================
# FEW PLOTS
# =============================================================================
x = np.linspace(7, 10.5)

plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], alpha=0.1)
plt.xlim(1.3, 2)
plt.ylim(1.3, 2)
plt.plot((0,10),(0,10),color='k')
plt.xlabel('redshift AD')
plt.ylabel('redshift BEAGLE')
plt.show()

plt.scatter(data['mass_AD_neb'], data['mass_BEAGLE_stellar'], alpha=0.1)
plt.xlim(7, 11)
plt.ylim(7, 11)
plt.plot((5,12),(5,12),color='k')
plt.xlabel('mass AD neb')
plt.ylabel('mass BEAGLE stellar')
plt.show()

plt.scatter(data['sfr_SAN'], data['sfr_BEAGLE_instant'], alpha=0.1)
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.plot((-3,3),(-3,3),color='k')
plt.xlabel('sfr SAN')
plt.ylabel('sfr BEAGLE instant')
plt.show()


# =============================================================================
# 2 sigma process
# =============================================================================

#color = data['msa_BEAGLE'] # colorbar for iterative 2 sigma clipping plots
#color = data['tau_BEAGLE']
#color = data['tauv_BEAGLE']
color = data['msa_BEAGLE'] - data['tau_BEAGLE']

if option_BEAGLE_mass == 0 and option_BEAGLE_sfr == 0:
    mass = data['mass_AD_neb'][data['sfr_SAN']>-90.0] # -99.0 is assigned when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)
    sfr = data['sfr_SAN'][data['sfr_SAN']>-90.0]
    color = color[data['sfr_SAN']>-90.0]
elif option_BEAGLE_mass == 1 and option_BEAGLE_sfr == 0:
    # SAME + mass_BEAGLE_stellar
    mass = data['mass_BEAGLE_stellar'][data['sfr_SAN']>-90.0]
    sfr = data['sfr_SAN'][data['sfr_SAN']>-90.0]
    color = color[data['sfr_SAN']>-90.0]
elif option_BEAGLE_mass == 0 and option_BEAGLE_sfr == 1:
    mass = data['mass_AD_neb']
    sfr = data['sfr_BEAGLE_instant']
elif option_BEAGLE_mass == 1 and option_BEAGLE_sfr == 1:   
    # SAME + sfr_BEAGLE_instant
    mass = data['mass_BEAGLE_stellar']
    sfr = data['sfr_BEAGLE_instant'] 
    


outliers = 1

# values saved as _0 prior to 2sigma clipping
mass0 = mass
sfr0 = sfr
fit0 = np.polyfit(mass0, sfr0, 1)
sfr_residuals0 = sfr0 - (fit0[0]*mass0 + fit0[1])
sigma0 = np.std(sfr_residuals0)

plt.scatter(mass0, sfr0)
plt.show()
print(len(mass))
    
while outliers > 0:
    
    fit = np.polyfit(mass, sfr, 1)
    sfr_residuals= sfr - (fit[0]*mass + fit[1])  
    sigma = np.std(sfr_residuals)
    idx = (abs(sfr_residuals)<2.0*sigma)    
    outliers = len(mass) - sum(idx)
    mass = mass[idx]
    sfr = sfr[idx]
    color = color[idx]
    
    # =============================================================================
    # PLOT ITERATIONS - colour coded
    # =============================================================================
    plt.scatter(mass, sfr, alpha=0.3, c=color)
    plt.xlim(7, 10)
    plt.ylim(-4, 3)
    plt.colorbar()
    plt.show()
    print(fit, sigma)   
    print('LEN MASS', len(mass))

    # =============================================================================
    # PLOT ITERATIONS - rising falling - color must be correctly set
    # =============================================================================
    idx_falling = (color > 0)
    idx_rising = (color < 0)
    plt.scatter(mass[idx_rising], sfr[idx_rising], alpha=0.3, color='k')
    plt.scatter(mass[idx_falling], sfr[idx_falling], alpha=0.3, color='r')
    plt.xlim(7, 10)
    plt.ylim(-4, 3)
    plt.show()
    print(fit, sigma)   
    print('LEN MASS', len(mass))
    
fit_clipping = fit


# =============================================================================
# FIRST YEAR REPORT
# =============================================================================

plt.figure(figsize=(6, 6))
plt.xlabel(r'$M_\mathrm{tot}$')
plt.ylabel(r'$\Psi$')
plt.scatter(mass0, sfr0, marker='x')
plt.scatter(mass, sfr, marker='x')
plt.plot(x, fit0[0]*x + fit0[1], color='#1f77b4', label='Without clipping', linewidth=2)
plt.plot(x, fit[0]*x + fit[1], color='#d62728', label='With clipping', linewidth=2)
plt.plot(x, beta_san[0]*x + alpha_san[0], color='#2ca02c', label='Santini+17', linewidth=2)
         
plt.plot(x, (beta_san[0]-beta_err_san[0])*x + (alpha_san[0]-alpha_err_san[0]), color='#2ca02c', linewidth=2, linestyle=':')
plt.plot(x, (beta_san[0]+beta_err_san[0])*x + (alpha_san[0]+alpha_err_san[0]), color='#2ca02c', linewidth=2, linestyle=':')
         

# report results
plt.plot(x, 0.84*x - 7.16, color='#9467bd', label='Report', linewidth=2) 
         
#plt.title(field.replace('_',''))
plt.xlim(7, 10.5)
plt.ylim(-2, 3)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
plt.show()

# =============================================================================
# FIRST YEAR REPORT - ZOOMED
# =============================================================================

plt.figure(figsize=(6, 6))
plt.xlabel(r'$M_\mathrm{tot}$')
plt.ylabel(r'$\Psi$')
#plt.scatter(mass0, sfr0, marker='x')
#plt.scatter(mass, sfr, marker='x')
plt.plot(x, fit0[0]*x + fit0[1], color='#1f77b4', label='Without clipping', linewidth=2)
plt.plot(x, fit[0]*x + fit[1], color='#d62728', label='With clipping', linewidth=2)
plt.plot(x, beta_san[0]*x + alpha_san[0], color='#2ca02c', label='Santini+17', linewidth=2)
         
plt.plot(x, (beta_san[0]-beta_err_san[0])*x + (alpha_san[0]-alpha_err_san[0]), color='#2ca02c', linewidth=2, linestyle=':')
plt.plot(x, (beta_san[0]+beta_err_san[0])*x + (alpha_san[0]+alpha_err_san[0]), color='#2ca02c', linewidth=2, linestyle=':')
         

# report results
plt.plot(x, 0.84*x - 7.16, color='#9467bd', label='Report', linewidth=2) 
         
#plt.title(field.replace('_',''))
plt.xlim(8.5, 9.0)
plt.ylim(-0.3, 0.4)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_santini.png')
plt.show()


# =============================================================================
# bootstrapping
# =============================================================================

mass = mass0
sfr = sfr0

iterations = 10 # 10000 for report
samples = len(mass)      # len(mass)== 1310, 1000 for report

alpha = []
beta = []
difference_8 = [] # difference between best fit and Santini at mass 10^8
difference_10 = []

for i in range(iterations):
    idx_random = np.random.choice(range(len(mass)), samples, replace=True)
#    print(np.sort(idx_random))
    mass_bs = mass[idx_random]
    sfr_bs = sfr[idx_random]
    
    outliers = 1
    
    while outliers > 0:
        
        fit = np.polyfit(mass_bs, sfr_bs, 1)
        sfr_residuals= sfr_bs - (fit[0]*mass_bs + fit[1])  
        sigma = np.std(sfr_residuals)
   
        idx = (abs(sfr_residuals)<2.0*sigma)    
        outliers = len(mass_bs) - sum(idx)
        mass_bs = mass_bs[idx]
        sfr_bs = sfr_bs[idx]
    
    alpha.append(fit[1])
    beta.append(fit[0])
    difference_8.append( (fit[1]+(8.0*fit[0])) - (alpha_san[0]+(8.0*beta_san[0])) )
    difference_10.append( (fit[1]+(10.0*fit[0])) - (alpha_san[0]+(10.0*beta_san[0])) )



# =============================================================================
# FIRST YEAR REPORT
# =============================================================================

#plt.title('alpha - intercept')
plt.xlabel(r'$\alpha$')
plt.ylabel('Count')
y, x, _ = plt.hist(alpha, bins=20)

plt.plot((np.median(alpha), np.median(alpha)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((np.percentile(alpha, 16), np.percentile(alpha, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
plt.plot((np.percentile(alpha, 84), np.percentile(alpha, 84)), (0, y.max()), color='#ff7f0e', linewidth=4, linestyle=':')
plt.plot((alpha_san[0], alpha_san[0]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
plt.plot((alpha_san[0]-alpha_err_san[0], alpha_san[0]-alpha_err_san[0]), (0, y.max()), label='68', color='#2ca02c', linewidth=4, linestyle=':')
plt.plot((alpha_san[0]+alpha_err_san[0], alpha_san[0]+alpha_err_san[0]), (0, y.max()), color='#2ca02c', linewidth=4, linestyle=':')
plt.plot((fit_clipping[1], fit_clipping[1]), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
plt.plot((-7.16, -7.16), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_alpha_bootstrap.png')
plt.show()

#plt.title('beta - slope')
plt.xlabel(r'$\beta$')
plt.ylabel('Count')
y, x, _ = plt.hist(beta, bins=20)

plt.plot((np.median(beta), np.median(beta)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((np.percentile(beta, 16), np.percentile(beta, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
plt.plot((np.percentile(beta, 84), np.percentile(beta, 84)), (0, y.max()), linewidth=4, linestyle=':')
plt.plot((beta_san[0], beta_san[0]), (0, y.max()), label='Santini+17', color='#2ca02c', linewidth=4)
plt.plot((beta_san[0]-beta_err_san[0], beta_san[0]-beta_err_san[0]), (0, y.max()), label='68', color='#2ca02c', linewidth=4, linestyle=':')
plt.plot((beta_san[0]+beta_err_san[0], beta_san[0]+beta_err_san[0]), (0, y.max()), color='#2ca02c', linewidth=4, linestyle=':')
plt.plot((fit_clipping[0], fit_clipping[0]), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
plt.plot((0.84, 0.84), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
plt.show()


#plt.title('difference_8')
plt.xlabel('difference8')
plt.ylabel('Count')
y, x, _ = plt.hist(difference_8, bins=20)

plt.plot((np.median(difference_8), np.median(difference_8)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((np.percentile(difference_8, 16), np.percentile(difference_8, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
plt.plot((np.percentile(difference_8, 84), np.percentile(difference_8, 84)), (0, y.max()), linewidth=4, linestyle=':')

plt.plot(((fit_clipping[1]+(8.0*fit_clipping[0])) - (alpha_san[0]+(8.0*beta_san[0])), (fit_clipping[1]+(8.0*fit_clipping[0])) - (alpha_san[0]+(8.0*beta_san[0]))), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
plt.plot(((-7.16+(8.0*0.84)) - (alpha_san[0]+(8.0*beta_san[0])), (-7.16+(8.0*0.84)) - (alpha_san[0]+(8.0*beta_san[0]))), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
plt.show()

#plt.title('difference_10')
plt.xlabel('difference10')
plt.ylabel('Count')
y, x, _ = plt.hist(difference_10, bins=20)

plt.plot((np.median(difference_10), np.median(difference_10)), (0, y.max()), label='Median', color='#ff7f0e', linewidth=4)
plt.plot((np.percentile(difference_10, 16), np.percentile(difference_10, 16)), (0, y.max()), label='68', color='#ff7f0e', linewidth=4, linestyle=':')
plt.plot((np.percentile(difference_10, 84), np.percentile(difference_10, 84)), (0, y.max()), linewidth=4, linestyle=':')

plt.plot(((fit_clipping[1]+(10.0*fit_clipping[0])) - (alpha_san[0]+(10.0*beta_san[0])), (fit_clipping[1]+(10.0*fit_clipping[0])) - (alpha_san[0]+(10.0*beta_san[0]))), (0, y.max()), label='Fit with clipping', color='k', linewidth=4)
plt.plot(((-7.16+(10.0*0.84)) - (alpha_san[0]+(10.0*beta_san[0])), (-7.16+(10.0*0.84)) - (alpha_san[0]+(10.0*beta_san[0]))), (0, y.max()), label='Report', color='#9467bd', linewidth=2)
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
plt.show()



print(len(mass))
















