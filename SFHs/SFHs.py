#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 17:57:05 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

plt.figure(figsize=(12,10))

t       = np.arange(0.1, 15, 0.1)
m       = 1

### 0 ### Iyer ### TopHat

A      = 1
t0     = 1
tau    = 6

sfr    = A * np.heaviside(t-t0, 1) * (1 - np.heaviside((t-t0-tau), 1))

# all a bit strange because discontinuous

m = 0           # t-t0 < 0
m = A*(t-t0)    #        0 < t-t0 < tau
m = A*tau       #                   tau < t-t0

A0 = 'anything'
A0 = m / (t-t0)
A0 = m / tau

plt.plot(t, sfr, label='TopHat')

### 1 ### Iyer ### Exponential

A      = 1
t0     = 2
tau    = 1

sfr    = A * np.heaviside(t-t0, 1) * np.exp(-(t-t0)/tau)
A0 = m / (tau * (1-(np.exp(-(t-t0)/tau)))) # t-t0 > 0

plt.plot(t, sfr, label='Exponential')

### 2 ### Carnall ### Delayed Exponential

A      = 1
t0     = 3
tau    = 1

sfr    = A * np.heaviside(t-t0, 1) * (t-t0) * np.exp(-(t-t0)/tau)
A0 = m / (-tau* ( ((np.exp(-(t-t0)/tau))*(t-t0+tau)) -tau) )


plt.plot(t, sfr, label='Delayed Exponential')

### 3 ### Santini ### Delayed Exponential

A      = 1
t0     = 4
tau    = 1

sfr    = A * np.heaviside(t-t0, 1) * ( ((t-t0)**2)/(tau**3) ) * np.exp(-(t-t0)/tau)
A0 = (-m*(tau**2)) / (((np.exp(-(t-t0)/tau))*((t0**2)-(2*t0*(tau+t))+(2*(tau**2))+(2*tau*t)+(t**2)))-(2*(tau**2)))

plt.plot(t, sfr, label='Santini Delayed Exponential')

### 4 ### Iyer ### Guassian

A      = 1
tpeak  = 5
tau    = 1

sfr    = A * np.exp( -(((t-tpeak)**2) / (2*(tau**2))) )
A0 = m / (-tau*((np.pi/2)**0.5)*((erf(-(t-tpeak)/((2**0.5)*tau)))-(erf(tpeak/((2**0.5)*tau)))))

plt.plot(t, sfr, label='Guassian')

### 5 ### Iyer ### Lognormal

A      = 5
t0     = 6
tau    = 1
t      = np.arange(t0+0.01, 15, 0.1)

sfr    = A * np.heaviside(t-t0, 1) * (1/t) * np.exp( -(((np.log(t-t0))**2) / (2*(tau**2))) )

plt.plot(t, sfr, label='Lognormal')



### 6 ### Carnall ### Double Power Law

A      = 1
alpha  = 10
beta   = 10
tau    = 10
t      = np.arange(0.1, 15, 0.1)

sfr    = A / (((t/tau)**alpha)+((t/tau)**-beta))

plt.plot(t, sfr, label='Double Power Law')

### PLOT ###

plt.xlabel('Age of Universe / Gyr')
plt.ylabel('SFR')

plt.legend()
plt.show()





### 6 ### Carnall ### Lognormal
'''
A      = 1
t0     = 7
tau    = 1
t      = np.arange(t0+0.01, 15, 0.1)

sfr    = A * (1/t) * np.exp( -(((np.log(t)-t0)**2) / (2*(tau**2))) )
plt.show()
plt.plot(t, sfr)
'''
### 8 ### Gladders ### Lognormal
'''
A      = 1
t0     = 8
tau    = 1
t      = np.arange(t0+0.01, 15, 0.1)

sfr    = A * (1/t) * np.exp( -(((np.log(t)-t0)**2) / (2*(tau**2))) )
plt.show()
plt.plot(t, sfr)
'''









































