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

revision = 'LE'
fsize=2
    
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}_combined.fits'.format(revision)
data_fits = fits.open(fileName)
#    print(data_fits[1].header)

ID = np.array(data_fits[1].data['id_1']).astype(float)

mass_in = 10**(data_fits[1].data['mass'])
msa_in = 10**data_fits[1].data['max_stellar_age']
tau_in = 10**(data_fits[1].data['tau'])
tau_exp_in = 10**(data_fits[1].data['tau_exp'])

mass_out = 10**(data_fits[1].data['mass_mean'])
msa_out = 10**data_fits[1].data['max_stellar_age_mean']
tau_out = 10**(data_fits[1].data['tau_mean'])
tau_exp_out = 10**(data_fits[1].data['tau_exp_mean'])

data_fits.close()

# calculated using python file 8
ID_rr_LE = np.array([ 6,  7,  8, 11, 15, 19, 20, 21, 25, 26, 27, 28, 31, 32, 34, 35, 36, 37, 38, 44, 46, 49, 50, 51, 52, 53, 55, 57, 58, 60, 63, 64, 66, 67, 68, 70, 71, 75, 76, 77, 79, 80, 81, 82, 83, 84, 85, 87, 88, 89, 94, 96, 97])
ID_ff_LE = np.array([22, 24, 43])
ID_rf_LE = np.array([69])
ID_fr_LE = np.array([  1,   2,   3,   4,   5,   9,  10,  12,  13,  14,  16,  17,  18,  23,  30,  33,  39,  40,  41,  42,  45,  48,  54,  56,  59,  61,  65,  73,  74,  78,  86,  90,  91,  93,  95,  98,  99, 100])

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s

xlin = np.linspace(1, 1e10, 100000)

for j in ID_ff_LE:
    
    i = (np.abs(ID - j)).argmin()
    
    plt.figure(figsize=(4*fsize, fsize))
    plt.title('{} ID {}'.format(revision.replace('_', ' '), j))
    plt.xlim(0, 1e10)
#    plt.xlim(0.9*msa, 1.1*msa)
    plt.ylim(0.0, 1.1)

    sfr_in = 1 * ((xlin-(ageUniv2-msa_in[i]))*np.heaviside(tau_in[i]-(xlin-(ageUniv2-msa_in[i])), 0) + tau_in[i]*np.exp((tau_in[i]-(xlin-(ageUniv2-msa_in[i])))/tau_exp_in[i])*np.heaviside((xlin-(ageUniv2-msa_in[i]))-tau_in[i], 1))
    plt.plot(xlin, sfr_in/max(sfr_in), label='INPUT SFH')
    
    sfr_out = 1 * ((xlin-(ageUniv2-msa_out[i]))*np.heaviside(tau_out[i]-(xlin-(ageUniv2-msa_out[i])), 0) + tau_out[i]*np.exp((tau_out[i]-(xlin-(ageUniv2-msa_out[i])))/tau_exp_out[i])*np.heaviside((xlin-(ageUniv2-msa_out[i]))-tau_out[i], 1))
    plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
    
    plt.plot((ageUniv2, ageUniv2), (0, 1))
    plt.legend()
    plt.show()
    
    print(mass_in[i], msa_in[i], tau_in[i], tau_exp_in[i])
    print(mass_out[i], msa_out[i], tau_out[i], tau_exp_out[i])

    


































