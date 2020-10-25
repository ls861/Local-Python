#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# Santini approx 1711 sources 1.3 < z < 6.0

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
sfr_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/calculate_santini_sfr_results/sfr_santini.npy'

AD = np.load(AD_location)
sfr_SAN = np.load(sfr_SAN_location)
print(AD.dtype.names)
print(len(AD))


# =============================================================================
# clusters only?
# =============================================================================
#idx = (AD['field']%2.0==0.0)
#AD = AD[idx]
#sfr_SAN = sfr_SAN[idx]
#print(len(AD))

# =============================================================================
# H band mag cut of < 27.5
# =============================================================================
idx = (AD['H160']<27.5)
AD = AD[idx]
sfr_SAN = sfr_SAN[idx]
print(len(AD))

# =============================================================================
# Photo-z from median of 6 methods
# =============================================================================
zLow = 1.3
zHigh = 2.0
wLim = (zHigh - zLow) / 2.0
zLim = zLow + wLim

idx = (abs(AD['ZBEST']-zLim) < wLim)
AD = AD[idx]
sfr_SAN = sfr_SAN[idx]
print(len(AD))

# =============================================================================
# Need sample to be complete above given mass (magnification NOT already included)
# =============================================================================
idx = ((np.log10(AD['MASTAR_NEB']*1e9)) > 8.4)
AD = AD[idx]
sfr_SAN = sfr_SAN[idx]
print(len(AD))

# =============================================================================
# Adding RELFLAG just to be sure...
# =============================================================================
idx = (AD['RELFLAG']==1.0)
AD = AD[idx]
sfr_SAN = sfr_SAN[idx]
print('RELFLAG==1', len(AD))

# =============================================================================
# Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
# =============================================================================
idx = (np.log10(AD['MASTAR_NEB']*1e9) - np.log10(AD['MAGNIF']) < 10.2)
AD = AD[idx]
sfr_SAN = sfr_SAN[idx]
print(len(AD))

# =============================================================================
# get some key columns
# =============================================================================
field_AD = AD['field']
id_AD = AD['ID']
mag_AD = AD['MAGNIF']
redshift_AD = AD['ZBEST']
mass_AD = np.log10(AD['MSTAR']*1e9) - np.log10(mag_AD)
mass_AD_neb = np.log10(AD['MASTAR_NEB']*1e9) - np.log10(mag_AD)
sfr_AD = np.log10(AD['SFR']) - np.log10(mag_AD)
sfr_AD_neb = np.log10(AD['SFR_NEB']) - np.log10(mag_AD)
sfr_SAN = sfr_SAN - np.log10(mag_AD)

# NO BEAGLE OUTPUTS USED YET
field_AD = np.array(field_AD)
id_AD = np.array(id_AD)
mag_AD = np.array(mag_AD)
redshift_AD = np.array(redshift_AD)
mass_AD = np.array(mass_AD)
mass_AD_neb = np.array(mass_AD_neb)
sfr_AD = np.array(sfr_AD)
sfr_AD_neb = np.array(sfr_AD_neb)
sfr_SAN = np.array(sfr_SAN)

# =============================================================================
# BEAGLE OUTPUTS
# =============================================================================
id_BEAGLE = []
mass_BEAGLE_tot = []
mass_BEAGLE_stellar = []
sfr_BEAGLE_instant = []
redshift_BEAGLE = []

'''
for i, field in enumerate(field_AD):

    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_inputs/BEAGLE_input_{}.fits'.format(int(field))
    data_fits = fits.open(fileName)
#    print(data_fits.info())
#    print(data_fits[1].header)
    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    data_fits.close()

#    GET BEAGLE ID FROM ORIGINAL

    if len(np.where(id_original==id_AD[i])[0]) == 0:
        print('AD OBJECT WAS NOT A BEAGLE INPUT')
        id_BEAGLE.append(-101)
    else:
        id_BEAGLE.append(id_input[np.where(id_original==id_AD[i])[0][0]])
        
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_summary_catalogues/BEAGLE_summary_catalogue_{}.fits'.format(int(field))
    data_fits = fits.open(fileName)
#    print(data_fits.info())
#    print(data_fits[2].header)
    id_input = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    
#    GET BEAGLE MASS, SFR and Z
    
    if len(np.where(id_input==id_BEAGLE[i])[0]) == 0:
        print('AD OBJECT WAS NOT FITTED BY BEAGLE')
        mass_BEAGLE_tot.append(-102.0)
        mass_BEAGLE_stellar.append(-102.0)
        redshift_BEAGLE.append(-102.0)     
        sfr_BEAGLE_instant.append(-102.0)
        data_fits.close()
    else:  
        idx = np.where(id_input==id_BEAGLE[i])[0][0]
        mass_BEAGLE_tot.append(data_fits['GALAXY PROPERTIES'].data['M_tot_median'][idx])
        mass_BEAGLE_stellar.append(data_fits['GALAXY PROPERTIES'].data['M_star_median'][idx])
        redshift_BEAGLE.append(data_fits['GALAXY PROPERTIES'].data['redshift_median'][idx])
        data_fits.close()
    
        data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_instant_SFRs/{}_instant_sfr_medians.p'.format(int(field)),'r'))
    #    print(data.keys())
        idx_sort = np.argsort(np.asarray(data['ID'], float))
        idx = np.where(np.asarray(data['ID'], float)[idx_sort]==id_BEAGLE[i])[0][0]
        sfr_BEAGLE_instant.append(data['log_instant_sfr_median.npy'][idx_sort][idx])
        
id_BEAGLE = np.array(id_BEAGLE)
mass_BEAGLE_tot = np.array(mass_BEAGLE_tot)
mass_BEAGLE_stellar = np.array(mass_BEAGLE_stellar)
sfr_BEAGLE_instant = np.array(sfr_BEAGLE_instant)
redshift_BEAGLE = np.array(redshift_BEAGLE)

np.save('id_BEAGLE',id_BEAGLE)
np.save('mass_BEAGLE_tot',mass_BEAGLE_tot)
np.save('mass_BEAGLE_stellar',mass_BEAGLE_stellar)
np.save('sfr_BEAGLE_instant',sfr_BEAGLE_instant)
np.save('redshift_BEAGLE',redshift_BEAGLE)

'''


id_BEAGLE = np.load('id_BEAGLE.npy')
mass_BEAGLE_tot = np.load('mass_BEAGLE_tot.npy')
mass_BEAGLE_stellar = np.load('mass_BEAGLE_stellar.npy')
sfr_BEAGLE_instant = np.load('sfr_BEAGLE_instant.npy')
redshift_BEAGLE = np.load('redshift_BEAGLE.npy')

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
# NOTES
# =============================================================================

x = np.linspace(7, 10.5)

'''
ALL THE FOLLOWING SHOULD BE EQUAL IN LENGTH:
    
field_AD
id_AD
mag_AD
redshift_AD
mass_AD
mass_AD_neb
sfr_AD
sfr_AD_neb
sfr_SAN

id_BEAGLE
mass_BEAGLE_tot
mass_BEAGLE_stellar
sfr_BEAGLE_instant
redshift_BEAGLE

NOTE sfr_SAN has -99.0 when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)

AND BEAGLE has -101.0 when the object was not a BEAGLE input and -102.0 if the object was a BEAGLE input, but was NOT fitted
'''

    
#%%

# =============================================================================
# 2 sigma process
# =============================================================================

mass = mass_AD_neb[sfr_SAN>-98.0] # -99.0 is assigned when not enough fluxes are +ve to calculate an SFR (or not enough filters in band, not the case for redshift > 1.3 though)
sfr = sfr_SAN[sfr_SAN>-98.0]

outliers = 1

mass0 = mass
sfr0 = sfr
fit0 = np.polyfit(mass0, sfr0, 1)
sfr_residuals0= sfr0 - (fit0[0]*mass0 + fit0[1])
sigma0 = np.std(sfr_residuals0)

plt.scatter(mass0, sfr0)
plt.show()
print(fit0, sigma0)   
print('LEN MASS', len(mass0))
    
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










