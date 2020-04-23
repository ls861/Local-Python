#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 18:09:05 2020

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits

param = 'DE'
revision = '108'
ID = 111


# =============================================================================
# OUTPUT - get BEAGLE parameters (<100)
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_100/fit_{}_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(revision, param)
data_fits = fits.open(fileName)

id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int) - 1

sfr_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_mean'])
sfr_68_b1 = np.log10(data_fits['STAR FORMATION'].data['SFR_68.00'])

ssfr_b1 = data_fits['STAR FORMATION'].data['sSFR_mean']
ssfr_68_b1 = data_fits['STAR FORMATION'].data['sSFR_68.00']

mass_b1 = data_fits['POSTERIOR PDF'].data['mass_mean']
msa_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_mean']
tauV_eff_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_mean']
metallicity_b1 = data_fits['POSTERIOR PDF'].data['metallicity_mean']
tau_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_mean']

mass_68_b1 = data_fits['POSTERIOR PDF'].data['mass_68.00']
msa_68_b1 = 10**data_fits['POSTERIOR PDF'].data['max_stellar_age_68.00']
tauV_eff_68_b1 = data_fits['POSTERIOR PDF'].data['tauv_eff_68.00']
metallicity_68_b1 = data_fits['POSTERIOR PDF'].data['metallicity_68.00']
tau_68_b1 = 10**data_fits['POSTERIOR PDF'].data['tau_68.00']

if revision in ['999']:
    nebular_logU_b1 = np.full(len(id_b1), -2.5)
    nebular_logU_68_b1 = np.full((len(id_b1),2), -2.5)
else:
    nebular_logU_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_mean']        
    nebular_logU_68_b1 = data_fits['POSTERIOR PDF'].data['nebular_logu_68.00']      

if revision in ['108']:
    nebular_xi_b1 = np.full(len(id_b1), 0.3)
    nebular_xi_68_b1 = np.full((len(id_b1),2), 0.3)
else:
    nebular_xi_b1 = data_fits['POSTERIOR PDF'].data['nebular_xi_mean']      
    nebular_xi_68_b1 = data_fits['POSTERIOR PDF'].data['nebular_xi_68.00']
    
redshift_b1 = data_fits['POSTERIOR PDF'].data['redshift_mean']
redshift_68_b1 = data_fits['POSTERIOR PDF'].data['redshift_68.00']

H160_b1 = (23.9 - data_fits['APPARENT MAGNITUDES'].data['HST_WFC3_IR_F160W_APP_mean']) / 2.5
#H160_68_b1 = data_fits['APPARENT MAGNITUDES'].data['HST_WFC3_IR_F160W_APP_68.00']

data_fits.close()

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/ASTRODEEP/{}_{}_chi2.fits'.format(revision, param)
data_fits = fits.open(fileName)
chi2 = data_fits[1].data['chi2']  
data_fits.close()


# =============================================================================
# INPUT
# =============================================================================

#header = ['b_B435', 'b_errB435', 'b_V606', 'b_errV606', 'b_I814', 'b_errI814', 'b_Y105', 'b_errY105', 'b_J125', 'b_errJ125', 'b_JH140', 'b_errJH140', 'b_H160', 'b_errH160', 'b_Ks', 'b_errKs', 'b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2', 'ZBEST', 'field', 'ID']

D1_new = np.load('D1_new.npy')

ZBEST = D1_new['ZBEST'][id_b1]
b_H160 = D1_new['b_H160'][id_b1]



# =============================================================================
# Timing Tests
# =============================================================================

runtime = np.empty(len(id_b1))

for i, ID in enumerate(id_b1 + 1):

    data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_{}/{}_BEAGLE.fits.gz'.format(revision, param, ID))
    #print(data_fits.info())

    runtime[i] = data_fits[0].header['HIERARCH RUN TIME'] / 3600
    data_fits.close()

    
# =============================================================================
# PLOTS
# =============================================================================

plt.title('runtime')
plt.hist(runtime)
plt.show()

plt.title('input Z vs runtime')
plt.scatter(ZBEST, runtime, marker='x', c=chi2, norm=colors.LogNorm())
plt.colorbar(label='chi2')
plt.show()

plt.title('output Z vs runtime')
plt.scatter(redshift_b1, runtime, marker='x', c=chi2, norm=colors.LogNorm())
plt.colorbar(label='chi2')
plt.show()

plt.title('in vs out Z')
plt.scatter(ZBEST, redshift_b1, marker='x', c=chi2, norm=colors.LogNorm())
plt.plot((0, 15), (0, 15))
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.colorbar(label='chi2')
plt.show()

plt.title('input H160 vs runtime')
plt.scatter(np.log10(abs(b_H160[b_H160>=0])+1e-6), runtime[b_H160>=0], marker='x')
plt.scatter(np.log10(abs(b_H160[b_H160<0])+1e-6), runtime[b_H160<0], color='r', marker='x')
plt.show()

plt.title('output H160 vs runtime')
plt.scatter(H160_b1, runtime, marker='x', c=chi2, norm=colors.LogNorm())
plt.colorbar(label='chi2')
plt.show()

plt.title('in vs out H160')
plt.scatter(np.log10(abs(b_H160)+1e-6), H160_b1, marker='x', c=chi2, norm=colors.LogNorm())
plt.plot((-9, 3), (-9, 3))
plt.xlim(-9, 3)
plt.ylim(-9, 3)
plt.colorbar(label='chi2')
plt.show()


plt.title('in Z vs in H160')
plt.scatter(ZBEST, np.log10(abs(b_H160)+1e-6), marker='x', c=runtime)
#plt.plot((-9, 3), (-9, 3))
#plt.xlim(-9, 3)
#plt.ylim(-9, 3)
plt.colorbar(label='runtime (hours)')
plt.show()













