#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:59:24 2020

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.special import erf
from scipy.integrate import quad
import matplotlib.colors as mcolors
import math


cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)
z           = np.random.uniform(low=4.5, high=5.5, size=size)                   # redshift

age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = 10**np.random.uniform(low=6, high=12, size=size)                  # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation
tau_arr     = 10**np.random.uniform(low=8.5, high=10, size=size)                   # yr, width of function

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

alpha_arr   = 10**np.random.uniform(low=-1, high=2, size=size)                  # cannot make upper limit 3 due to beta_max
beta_arr    = 10**np.random.uniform(low=-1, high=2, size=size)



k = 100                                                                         # the factor in beta_min
#beta_max    = -(np.log10((2*k) - ((t_arr/tau_arr)**(alpha_arr)))) / np.log10(t_arr/tau_arr) # limits left slope

'''
beta_max   = np.zeros(size)

for i in range(size):
    print(t_arr[i], tau_arr[i], alpha_arr[i])
    beta_max[i]    = -(np.log10((2*k) - ((t_arr[i]/tau_arr[i])**(alpha_arr[i])))) / np.log10(t_arr[i]/tau_arr[i]) # limits left slope
'''

#test = 0
'''
for i in range(size):
    if math.isnan(beta_max[i]):
        test+=1

    
    if beta_max[i] > 100 or math.isnan(beta_max[i]):
        beta_max[i] = np.log10(100)
        
    elif beta_max[i] < 0.1:
        beta_max[i] = np.log10(0.1)
        
    else:
        beta_max[i] = np.log10(beta_max[i])

print('nan', test)

beta_arr    = 10**np.random.uniform(low=-1, high=beta_max, size=size)
'''




### Delayed Exponential ###

sfr_arr = []
mass_used = []

ssfr_arr = []
tau_log10_arr = []

for i in range(size):

    if age_z_15 < t0_arr[i]:
    
        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0       = t0_arr[i]                                                    # yrs        
        tau     = tau_arr[i]                                                    # yrs
        alpha   = alpha_arr[i]
        beta    = beta_arr[i]
        
        integrand = lambda T: 1 / (((T/tau)**alpha)+((T/tau)**-beta))
        integral  = quad(integrand, 0, t)
        
        A = m / integral[0] 
        sfr = A / (((t/tau)**alpha)+((t/tau)**-beta))
        
        mass_used.append(m_arr[i])
        sfr_arr.append(np.log10(sfr))
        
        ssfr_arr.append(np.log10(sfr/m))
        tau_log10_arr.append(np.log10(tau))

        
        xax = (10**9)*np.arange(0.1, 15, 0.01)
        upper = max(A / (((xax/tau)**alpha)+((xax/tau)**-beta)))
        '''
        plt.figure(figsize=(10, 5))
        plt.plot(xax, A / (((xax/tau)**alpha)+((xax/tau)**-beta)))
        plt.plot((t, t), (0, upper))
        #plt.plot((t0, t0), (0, upper))
        plt.xlim((0, 14E9))
        plt.ylim((0, upper))
        plt.show()
        print(A, m, t, t0, alpha, beta)
        '''
            
print(len(t_arr))
print(len(sfr_arr))

### ### ### ### ###

plt.figure(figsize=(7, 5))
plt.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")

#plt.xlim(6,11)
#plt.ylim(-4,3)

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s

plt.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(1E6)))

plt.show()

### ### ### ### ###

plt.figure(figsize=(7, 5))
plt.hist2d(ssfr_arr, tau_log10_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log ssfr")
plt.ylabel("log tau")

#plt.xlim(6,11)
#plt.ylim(-4,3)

plt.show()

### ### ### ### ###