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

# Santini approx 1711 sources 1.3 < z < 6.0

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
                'mag_AD':               AD['MAGNIF'], 
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
    
                }

# =============================================================================
# clusters only?
# =============================================================================
idx1 = (AD['field']%1.0==0.0) # all
#idx1 = (AD['field']%2.0==0.0) # clusters

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
idx4 = ((np.log10(AD['MASTAR_NEB']*1e9)) > 8.4)
idx = np.logical_and(idx,idx4)

# =============================================================================
# Adding RELFLAG just to be sure...
# =============================================================================
idx5 = (AD['RELFLAG']==1.0)
idx = np.logical_and(idx,idx5)

# =============================================================================
# Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
# =============================================================================
idx6 = (np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']) < 10.2)
idx = np.logical_and(idx,idx6)

# =============================================================================
# COMBINE THE ABOVE
# =============================================================================

for key in data.keys():
    print(key, len(data[key]), min(data[key]))
#    print(data[key])

for key in data.keys():
    data[key] = data[key][idx]
    print(key, len(data[key]), min(data[key]))

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

plt.scatter(data['redshift_AD'], data['redshift_BEAGLE'], alpha=0.01)
plt.xlim(0, 10)
plt.ylim(0, 10)
plt.plot((0,10),(0,10),color='k')
plt.show()

plt.scatter(data['mass_AD_neb'], data['mass_BEAGLE_stellar'], alpha=0.01)
plt.xlim(5, 12)
plt.ylim(5, 12)
plt.plot((5,12),(5,12),color='k')
plt.show()

plt.scatter(data['sfr_AD_neb'], data['sfr_BEAGLE_instant'], alpha=0.01)
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.plot((-3,3),(-3,3),color='k')
plt.show()


#%%

# =============================================================================
# 2 sigma process
# =============================================================================

mass = data['mass_AD_neb'][data['sfr_SAN']>-90.0] # -99.0 is assigned when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)
sfr = data['sfr_SAN'][data['sfr_SAN']>-90.0]

outliers = 1

mass0 = mass
sfr0 = sfr
fit0 = np.polyfit(mass0, sfr0, 1)
sfr_residuals0 = sfr0 - (fit0[0]*mass0 + fit0[1])
sigma0 = np.std(sfr_residuals0)

plt.scatter(mass0, sfr0)
plt.show()
    
while outliers > 0:
    
    fit = np.polyfit(mass, sfr, 1)
    sfr_residuals= sfr - (fit[0]*mass + fit[1])  
    sigma = np.std(sfr_residuals)
    idx = (abs(sfr_residuals)<2.0*sigma)    
    outliers = len(mass) - sum(idx)
    mass = mass[idx]
    sfr = sfr[idx]    
    plt.scatter(mass, sfr)
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









'''
# =============================================================================
# bootstrapping
# =============================================================================
#%%


mass = mass_AD_neb[sfr_SAN>-98.0] # -99.0 is assigned when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)
sfr = sfr_SAN[sfr_SAN>-98.0]


print(len(mass_AD_neb))

iterations = 1000 # 10000 for report
samples = 100      # len(mass)== 1310, 1000 for report

alpha = []
beta = []

for i in range(iterations):
    idx_random = np.random.choice(range(len(mass)), samples, replace=True)
#    print(np.sort(idx_random))
    mass_bs = mass[idx_random]
    sfr_bs = sfr[idx_random]
    
    outliers = 1
    
    mass1 = mass_bs
    sfr1 = sfr_bs
    fit1 = np.polyfit(mass1, sfr1, 1)
    sfr_residuals1= sfr1 - (fit1[0]*mass1 + fit1[1])
    sigma1 = np.std(sfr_residuals1)
    
    while outliers > 0:
        
        fit = np.polyfit(mass_bs, sfr_bs, 1)
        sfr_residuals= sfr_bs - (fit[0]*mass_bs + fit[1])  
        sigma = np.std(sfr_residuals)
   
        idx = (abs(sfr_residuals)<2.0*sigma)    
        outliers = len(mass_bs) - sum(idx)
        mass_bs = mass_bs[idx]
        sfr_bs = sfr_bs[idx]    

    # logSFR = alpha log(M / M_9.7) + beta
    z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
    A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
    A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
    B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
    B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))
    
    alpha_san = B_san - 9.7*A_san
    beta_san = A_san
    
    alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
    beta_err_san = A_err_san
    
    x = np.linspace(min(mass_bs), max(mass_bs))
    
    alpha.append(fit[1])
    beta.append(fit[0])



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
plt.legend()
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/334_beta_bootstrap.png')
plt.show()

'''








