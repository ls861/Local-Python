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
from scipy.integrate import quad
import math

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E6)
z           = np.random.uniform(low=4.5, high=5.5, size=size)                   # redshift

age_z_4p5   = cd.age(4.5, **cosmo)/cc.yr_s                                      # yr, (log, 9.107)
age_z_5p5   = cd.age(5.5, **cosmo)/cc.yr_s                                      # yr,
age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = 10**np.random.uniform(low=6, high=10, size=size)                  # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

### ### ### ### ### ### ###
### Delayed Exponential ###
### ### ### ### ### ### ### 

tau_arr     = 10**np.random.uniform(low=7.5, high=10, size=size)                # yr, width of function

sfr_arr = []
mass_used = []

ssfr_arr = []
tau_log10_arr = []

for i in range(size):

    if age_z_15 < t0_arr[i]: # this ensures galaxy formation begins after z=15.
    
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


print(len(t_arr))
print(len(sfr_arr))
 
### ### ### ### ### sfr vs mass

fig, ax = plt.subplots(figsize=(8, 4))

h = ax.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
fig.colorbar(h[3], ax=ax, label='weighting in prior')

ax.set_xlabel("log mass")
ax.set_ylabel("log sfr")

plt.plot(massArr, np.log10(10**massArr/(1E6)), label='t=1E6')
plt.plot(massArr, np.log10(10**massArr/age_z_4p5), label='z=4.5')
plt.plot(massArr, np.log10(10**massArr/age_z_5p5), label='z=5.5')
plt.legend()
plt.show()

### ### ### ### ###

 
### ### ### ### ### tau vs ssfr

fig, ax = plt.subplots(figsize=(8, 4))
h = ax.hist2d(ssfr_arr, tau_log10_arr, bins=[50,100], cmap='Blues', normed=True)
fig.colorbar(h[3], ax=ax, label='weighting in prior')

ax.set_xlabel("log ssfr")
ax.set_ylabel("log tau")
plt.show()

### ### ### ### ###

 

### ### ### ### ###   ###
### Double Power Law1 ###
### ### ### ### ###   ###       

tau_arr     = np.random.uniform(low=7E8, high=t_arr, size=size)                 # yr

alpha_arr   = (10**np.random.uniform(low=-1, high=3, size=size))
beta_arr    = (10**np.random.uniform(low=-1, high=3, size=size))

### ### ### ### ###

sfr_arr     = []
mass_used   = []
int_err     = []

print(age_z_15)
for i in range(size):

    if age_z_15 < t0_arr[i]: # age of galaxy is less than time available for star formation

        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0      = t0_arr[i]                                                     # yrs
        tau     = tau_arr[i]                                                    # yrs
        
        alpha   = alpha_arr[i]
        beta    = beta_arr[i]
        
        integrand = lambda T: 1 / (((T/tau)**alpha)+((T/tau)**-beta))
        integral  = quad(integrand, t0, t)
            
        if integral[0] > 0:                                                     # not zero
            A = m / integral[0]                                                 # solar masses / yr
            
            if A > 1E300:
                #pass
                print('lt', A, m, integral[0], (((t/tau)**alpha)+((t/tau)**-beta)), sfr)
            else:
                int_err.append(integral[1])
                sfr = A / (((t/tau)**alpha)+((t/tau)**-beta))
                
                if sfr < 1E-10:
                    pass
                    #print(('sfr', A, m, integral[0], (((t/tau)**alpha)+((t/tau)**-beta)), sfr))
                else:
                    mass_used.append(m_arr[i])                                      # log, mass of galaxy
                    sfr_arr.append(np.log10(sfr)) 
        else:
            print('TEST', integral[0])
            T = np.linspace(0, 1E10, 1000)
            plt.plot(T, 1 / (((T/tau)**alpha)+((T/tau)**-beta)))
            plt.show()


print(len(t_arr), len(sfr_arr), max(int_err))
print(max(mass_used), min(mass_used), max(sfr_arr), min(sfr_arr))

plt.figure(figsize=(8, 4))
plt.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")

#plt.xlim(6,11)
#plt.ylim(-4,3)

plt.plot(massArr, np.log10(10**massArr/(1E6)), label='t=1E6')
plt.plot(massArr, np.log10(10**massArr/age_z_4p5), label='z=4.5')
plt.plot(massArr, np.log10(10**massArr/age_z_5p5), label='z=5.5')
plt.legend()
plt.show()


### ### ### ### ###   ###
### Double Power Law2 ###
### ### ### ### ###   ###  


tau_arr     = 10**np.random.uniform(low=9, high=10, size=size)                  # yr, width of function

alpha_arr   = 10**np.random.uniform(low=-1, high=2, size=size)                  # cannot make upper limit 3 due to beta_max

k = 100                                                                         # the factor in beta_min

beta_max    = -(np.log10((2*k) - ((t_arr/tau_arr)**(alpha_arr)))) / np.log10(t_arr/tau_arr) # limits left slope



test = 0
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

### Delayed Exponential ###

sfr_arr = []
mass_used = []

ssfr_arr = []
tau_log10_arr = []
alpha_log10_arr = []
beta_log10_arr = []

for i in range(size):

    if age_z_15 < t0_arr[i]:
    
        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0       = 0                                                    # yrs        
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
        alpha_log10_arr.append(np.log10(alpha))
        beta_log10_arr.append(np.log10(beta))
        
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

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s

plt.plot(massArr, np.log10(10**massArr/(1E6)), label='t=1E6')
plt.plot(massArr, np.log10(10**massArr/ageUniv_4p5), label='z=4.5')
plt.plot(massArr, np.log10(10**massArr/ageUniv_5p5), label='z=5.5')

plt.show()

### ### ### ### ### TAU

plt.figure(figsize=(7, 5))
plt.hist2d(ssfr_arr, tau_log10_arr, bins=[np.arange(-10.5, -7.6, 0.05), np.arange(9, 10.1, 0.05)], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log ssfr")
plt.ylabel("log tau")

plt.show()

### ### ### ### ### ALPHA

plt.figure(figsize=(7, 5))
plt.hist2d(ssfr_arr, alpha_log10_arr, bins=[500, 100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log ssfr")
plt.ylabel("log alpha")

plt.show()

### ### ### ### ### BETA

plt.figure(figsize=(7, 5))
plt.hist2d(ssfr_arr, beta_log10_arr, bins=[500, 100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log ssfr")
plt.ylabel("log beta")

plt.show()

### ### ### ### ###

plt.figure(figsize=(7, 5))
plt.hist2d(tau_log10_arr, mass_used, bins=[50, 100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log tau")
plt.ylabel("log mass")

#plt.xlim(-10.5,-7.5)
#plt.ylim(9,10)

plt.show()

### ### ### ### ###










