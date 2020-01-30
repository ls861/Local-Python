#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:41:48 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)
z           = np.random.uniform(low=4.5, high=5.5, size=size)                   # redshift

age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = 10**np.random.uniform(low=6, high=10, size=size)                  # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation
tau_arr     = 10**np.random.uniform(low=7.5, high=10, size=size)                # yr, width of function

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)


### Delayed Exponential ###

sfr_arr = []
mass_used = []

ssfr_arr = []
tau_log10_arr = []

for i in range(size):

    if age_z_15 < t0_arr[i]:
    
        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0      = t0_arr[i]                                                     # yrs
        tau     = tau_arr[i]                                                    # yrs

        A = m / (-tau* ( ((np.exp(-(t-t0)/tau))*(-t0+tau+t)) - tau) )
        # https://www.wolframalpha.com/input/?i=%28t-B%29*exp%28-%28t-B%29%2FC%29  
        
        sfr = A * (t-t0) * np.exp(-(t-t0)/tau)
        mass_used.append(m_arr[i])
        sfr_arr.append(np.log10(sfr))
        
        ssfr_arr.append(np.log10(sfr/m))
        tau_log10_arr.append(np.log10(tau))

        '''
        xax = (10**9)*np.arange(0.1, 15, 0.01)
        plt.figure(figsize=(10, 5))
        plt.plot(xax, A * (xax-t0) * np.exp(-(xax-t0)/tau))
        plt.plot((t, t), (0, 1E100))
        plt.plot((t0, t0), (0, 1E10))
        plt.xlim((0, 14E9))
        plt.ylim((0, max(A * (xax-t0) * np.exp(-(xax-t0)/tau))))
        plt.show()
        print(A, m)
        '''
            
print(len(t_arr))
print(len(sfr_arr))

### ### ### ### ### using fig, ax instead of plt

fig, ax = plt.subplots(figsize=(10, 8))

h = ax.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
fig.colorbar(h[3], ax=ax, label='weighting in prior')

ax.set_xlabel("log mass")
ax.set_ylabel("log sfr")

ax.set_xlim(6,11)
ax.set_ylim(-4,3)

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s

ax.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
ax.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))
ax.plot(massArr, np.log10(10**massArr/(1E6)))

fig.show()


































