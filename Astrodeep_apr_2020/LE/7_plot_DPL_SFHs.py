#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 14:53:16 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_apr_2020/LE/mock_MS_parameters_100_LE.fits'
data_fits = fits.open(fileName)
#print(data_fits[1].header)
msa = 10**data_fits[1].data['max_stellar_age']
tau = 10**(data_fits[1].data['tau'])
tau_exp = 10**(data_fits[1].data['tau_exp'])
A = data_fits[1].data['A']
sfr = data_fits[1].data['sfr']
data_fits.close()

ageUniv2 = 3228839870.9122815
xlin = np.linspace(1, 1e10, 100000)
grad = np.empty(len(A))



plt.figure(figsize=(10, 10))
plt.xlim(0, 1e10)
plt.ylim(0, 1)

for i in range(len(A)):
#for i in [6]:
    
    # nice trick to find index in xlin which has value closest to ageUniv2
    idx = (np.abs(xlin-ageUniv2).argmin())
    
    # replace msa with (xlin-(ageUniv2-msa[i]))

    sfr_calc = A[i] * ((xlin-(ageUniv2-msa[i]))*np.heaviside(tau[i]-(xlin-(ageUniv2-msa[i])), 0) + tau[i]*np.exp((tau[i]-(xlin-(ageUniv2-msa[i])))/tau_exp[i])*np.heaviside((xlin-(ageUniv2-msa[i]))-tau[i], 1))
#    print(alpha[i], beta[i], tau[i], A[i], sfr[i], sfr_calc)
    plt.plot(xlin, sfr_calc/max(sfr_calc))
#    plt.plot(xlin, np.gradient(sfr_calc, xlin))
    grad[i] = np.gradient(sfr_calc, xlin)[idx]
    
plt.show()

plt.hist(grad, bins=50)
plt.show()


plt.hist(tau, bins=50)
plt.show()



plt.figure(figsize=(10, 10))
plt.xlim(0, 1e10)
plt.ylim(0.9, 1)

rising = np.array(range(len(A)))[grad<0]
print(rising)

for i in rising:
#for i in [6]:
    
    # nice trick to find index in xlin which has value closest to ageUniv2
    idx = (np.abs(xlin-ageUniv2).argmin())

    sfr_calc = A[i] * ((xlin-(ageUniv2-msa[i]))*np.heaviside(tau[i]-(xlin-(ageUniv2-msa[i])), 0) + tau[i]*np.exp((tau[i]-(xlin-(ageUniv2-msa[i])))/tau_exp[i])*np.heaviside((xlin-(ageUniv2-msa[i]))-tau[i], 1))
#    print(alpha[i], beta[i], tau[i], A[i], sfr[i], sfr_calc)
    plt.plot(xlin, sfr_calc/max(sfr_calc))
#    plt.plot(xlin, np.gradient(sfr_calc, xlin))
    grad[i] = np.gradient(sfr_calc, xlin)[idx]
    
plt.show()


plt.figure(figsize=(10, 10))
plt.xlim(0, 1e10)
plt.ylim(0.9, 1)


falling = np.array(range(len(A)))[grad>0]
print(falling)

for i in falling:
#for i in [6]:
    
    # nice trick to find index in xlin which has value closest to ageUniv2
    idx = (np.abs(xlin-ageUniv2).argmin())

    sfr_calc = A[i] * ((xlin-(ageUniv2-msa[i]))*np.heaviside(tau[i]-(xlin-(ageUniv2-msa[i])), 0) + tau[i]*np.exp((tau[i]-(xlin-(ageUniv2-msa[i])))/tau_exp[i])*np.heaviside((xlin-(ageUniv2-msa[i]))-tau[i], 1))
#    print(alpha[i], beta[i], tau[i], A[i], sfr[i], sfr_calc)
    plt.plot(xlin, sfr_calc/max(sfr_calc))
#    plt.plot(xlin, np.gradient(sfr_calc, xlin))
    grad[i] = np.gradient(sfr_calc, xlin)[idx]
    
plt.show()


