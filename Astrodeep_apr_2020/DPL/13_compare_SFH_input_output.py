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

revision = '012_010'
fsize=2
    
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/{}_input_and_output.fits'.format(revision)
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
ID_rr_012_010 = np.array([2, 40, 72])
ID_ff_012_010 = np.array([1, 3, 4, 6, 8, 12, 16, 18, 21, 22, 24, 25, 26, 27, 28, 29, 35, 41, 43, 44, 47, 49, 50, 51, 52, 53, 55, 57, 58, 61, 63, 66, 67, 68, 71, 73, 74, 76, 77, 79, 80, 82, 83, 84, 85, 90, 92, 93, 94, 95, 96, 98, 99])
ID_rf_012_010 = np.array([5, 7, 9, 10, 11, 13, 14, 32, 33, 34, 37, 38, 56, 62, 69, 70, 75, 78, 86, 87, 88, 89, 91, 97])
ID_fr_012_010 = np.array([17, 81, 100])

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s
msa = ageUniv2 - ageUniv999

xlin = np.linspace(1, 1e10, 100000)

for j in ID_rf_012_010:
    
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

    


































