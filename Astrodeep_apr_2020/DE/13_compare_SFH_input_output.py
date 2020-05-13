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

revision = 'DE'
fsize=2
    
fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/{}_combined.fits'.format(revision)
data_fits = fits.open(fileName)
#    print(data_fits[1].header)

ID = np.array(data_fits[1].data['id_1']).astype(float)

mass_in = 10**(data_fits[1].data['mass'])
msa_in = 10**data_fits[1].data['max_stellar_age']
tau_in = 10**(data_fits[1].data['tau'])

mass_out = 10**(data_fits[1].data['mass_mean'])
msa_out = 10**data_fits[1].data['max_stellar_age_mean']
tau_out = 10**(data_fits[1].data['tau_mean'])

data_fits.close()

# calculated using python file 8
ID_rr_DE = np.array([ 4,  5,  8, 10, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 34, 35, 36, 37, 39, 41, 42, 43, 44, 45, 47, 49, 50, 51, 53, 55, 56, 58, 61, 62, 63, 64, 65, 66, 68, 70, 71, 73, 75, 78, 79, 81, 83, 84, 88, 89, 90, 92, 93, 94, 95, 96, 98, 99])
ID_ff_DE = np.array([ 3, 59])
ID_rf_DE = np.array([])
ID_fr_DE = np.array([2,   7,   9,  15,  18, 33,  38,  40,  52,  69,  72,  74,  76,  77,  80,  82,  85,  86,87,  91,  97, 100])
test = [(4)]


cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
ageUniv2 = cd.age(2.0, **cosmo)/cc.yr_s
ageUniv999 = cd.age(999.0, **cosmo)/cc.yr_s

xlin = np.linspace(1, 1e10, 100000)

for j in test:
    
    i = (np.abs(ID - j)).argmin()
    
    plt.figure(figsize=(4*fsize, fsize))
    plt.title('{} ID {}'.format(revision.replace('_', ' '), j))
    plt.xlim(0, 1e10)
#    plt.xlim(0.9*msa, 1.1*msa)
    plt.ylim(0.0, 1.1)

    sfr_in = 1 * (xlin-(ageUniv2-msa_in[i]))*np.exp(-(xlin-(ageUniv2-msa_in[i]))/tau_in[i])
    plt.plot(xlin, sfr_in/max(sfr_in), label='INPUT SFH')
    
    sfr_out = 1 * (xlin-(ageUniv2-msa_out[i]))*np.exp(-(xlin-(ageUniv2-msa_out[i]))/tau_out[i])
    plt.plot(xlin, sfr_out/max(sfr_out), label='OUTPUT SFH')
    
    plt.plot((ageUniv2, ageUniv2), (0, 1))
    plt.legend()
    plt.show()
    
    print(mass_in[i], msa_in[i], tau_in[i])
    print(mass_out[i], msa_out[i], tau_out[i])

    


































