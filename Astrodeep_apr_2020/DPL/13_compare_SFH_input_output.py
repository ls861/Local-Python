#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 09:09:20 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cosmolopy.distance as cd
import cosmolopy.constants as cc

revision = 'DPL'
fsize=2
    
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}_combined.fits'.format(revision)
data_fits = fits.open(fileName)
#    print(data_fits[1].header)

ID = np.array(data_fits[1].data['id_1']).astype(float)

mass_in = 10**(data_fits[1].data['mass'])
alpha_in = data_fits[1].data['dpl_alpha']
beta_in = data_fits[1].data['dpl_beta']
tau_in = 10**(data_fits[1].data['tau'])

mass_out = 10**(data_fits[1].data['mass_mean'])
alpha_out = data_fits[1].data['dpl_alpha_mean']
beta_out = data_fits[1].data['dpl_beta_mean']
tau_out = 10**(data_fits[1].data['tau_mean'])

data_fits.close()

# calculated using python file 8
ID_rr_DPL = np.array([2])
ID_ff_DPL = np.array([1, 4, 6, 7, 10, 14, 15, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 41, 42, 43, 44, 47, 48, 50, 51, 52, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 74, 82, 83, 84, 86, 87, 88, 89, 90, 91, 92, 94, 97, 99, 100])
ID_rf_DPL = np.array([3, 5, 8, 11, 12, 13, 16, 23, 31, 39, 45, 46, 49, 53, 60, 64, 72, 73, 75, 76, 77, 79, 80, 81, 95, 96, 98])
ID_fr_DPL = np.array([33, 37, 93])

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s
msa = ageUniv2 - ageUniv999

xlin = np.linspace(1, 1e10, 100000)

for j in ID_rf_DPL:
    
    i = (np.abs(ID - j)).argmin()
    
    plt.figure(figsize=(4*fsize, fsize))
    plt.title('{} ID {}'.format(revision.replace('_', ' '), j))
    plt.xlim(0, 1e10)
#    plt.xlim(0.9*msa, 1.1*msa)
    plt.ylim(0.0, 1.1)

    sfr_in = 1 / (((xlin/tau_in[i])**alpha_in[i])+((xlin/tau_in[i])**-beta_in[i]))
    plt.plot(xlin, sfr_in/max(sfr_in), label='INPUT SFH')
    
    sfr_out = 1 / (((xlin/tau_out[i])**alpha_out[i])+((xlin/tau_out[i])**-beta_out[i]))
    plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
    
    plt.plot((msa, msa), (0, 1))
    plt.legend()
    plt.show()
    
    print(mass_in[i], alpha_in[i], beta_in[i], tau_in[i])
    print(mass_out[i], alpha_out[i], beta_out[i], tau_out[i])

    


































