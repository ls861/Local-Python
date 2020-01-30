#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:06:01 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

### ### MASS dependence on tau ### ###

### ### ### ### ### ### ###
### Delayed Exponential ###
### ### ### ### ### ### ### 

A = 1
t = 1.12 * 1E9              # time of observation / age of universe at z=5
msa = 5E8                   # 10**8.70
t0 = t - msa

tau_arr = 10**np.arange(7.5, 10.1, 0.1)
mass_arr = A*(-tau_arr * ( ((np.exp(-(t-t0)/tau_arr))*(-t0+tau_arr+t)) - tau_arr) )

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tau_arr, mass_arr / max(mass_arr))
ax.plot((msa, msa), (0, 1.1), color='k', linestyle=':')   
ax.set_xlim(1E10-t0, 0-t0)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TAU')
ax.set_ylabel('MASS / A')
plt.show()
    


### ### ### ### ###  ###
### Double Power Law ###
### ### ### ### ###  ###  

A = 300000000
t = 1.12 * 1E9              # time of observation / age of universe at z=5 / 10**9.049218
alpha = 3.
beta = 1.

tau_arr = 10**np.arange(7.5, 10.1, 0.1)
mass_arr = np.zeros(len(tau_arr))

tau_arr = 133000000.
t = 1E8
tau = 133000000.

'''
for i in range(len(tau_arr)):
    
    integrand = lambda T: 1 / (((T/tau_arr[i])**alpha)+((T/tau_arr[i])**-beta))
    integral  = quad(integrand, 0, t)


    mass_arr[i] = A * integral[0] 
'''


# https://www.wolframalpha.com/input/?i=%28%28t%2FC%29%5EB+%28B+-+A+%28t%2FC%29%5E%28A+%2B+B%29%29%29%2F%28t+%281+%2B+%28t%2FC%29%5E%28A+%2B+B%29%29%5E2%29&assumption=%22ClashPrefs%22+-%3E+%7B%22Math%22%7D

t_peak = np.exp( (   (alpha*np.log(tau)) - np.log(alpha) + (beta*np.log(tau)) + np.log(beta))    /  (alpha + beta))
tau_peak = np.exp(    ( (alpha*np.log(t)) + np.log(alpha) + (beta*np.log(t)) - np.log(beta))  /  (alpha + beta) )



'''
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(tau_arr, mass_arr / max(mass_arr))
ax.plot((t, t), (0, 1.1), color='k', linestyle=':')   
ax.plot((tau_peak, tau_peak), (0, 1.1), color='r', linestyle=':')   
ax.set_xlim(1E10, 0)
ax.set_ylim(0, 1.1)
ax.set_xlabel('TAU')
ax.set_ylabel('MASS / A')
plt.show()


print(tau_peak)
'''












