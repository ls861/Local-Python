#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:53:16 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_004.fits'
data_fits = fits.open(fileName)
#print(data_fits[1].header)
alpha = data_fits[1].data['dpl_alpha']
beta = data_fits[1].data['dpl_beta']
tau = 10**(data_fits[1].data['tau'])
A = data_fits[1].data['A']
sfr = data_fits[1].data['sfr']
data_fits.close()

ageUniv2 = 3228839870.9122815
xlin = np.linspace(1, 1e10, 100000)
grad = np.empty(len(A))

# nice trick to find index in xlin which has value closest to ageUniv2
idx = (np.abs(xlin - ageUniv2)).argmin()



plt.figure(figsize=(10, 10))
plt.xlim(0, 1e10)
#plt.ylim(0, 50)

for i in range(len(A)):
#for i in [0]:
    sfr_calc = A[i] / (((xlin/tau[i])**alpha[i])+((xlin/tau[i])**-beta[i]))
#    print(alpha[i], beta[i], tau[i], A[i], sfr[i], sfr_calc)
    plt.plot(xlin, sfr_calc/max(sfr_calc))
#    plt.plot(xlin, np.gradient(sfr_calc, xlin))
    grad[i] = np.gradient(sfr_calc, xlin)[idx]
    
plt.show()

plt.hist(grad, bins=50)
plt.show()


plt.hist(tau, bins=50)
plt.show()












