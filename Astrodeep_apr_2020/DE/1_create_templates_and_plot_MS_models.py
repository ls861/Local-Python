#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:09:12 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import truncnorm
import cosmolopy.distance as cd
import cosmolopy.constants as cc

nObj = 10000

intercept = -8.
slope = 1.
scatter = 0.3

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
zObs=2.0
ageUniv = cd.age(zObs, **cosmo)/cc.yr_s
ageUniv5 = cd.age(5.0, **cosmo)/cc.yr_s

msa = 10**np.random.uniform(low=6, high=np.log10(ageUniv), size=nObj)                      # yrs
tau = 10**np.random.uniform(low=7, high=10.5, size=nObj)                      # yrs
mass = 10**np.random.uniform(low=7, high=12, size=nObj)                         # solar masses



A = np.zeros(nObj)

for i in range(nObj):
    integrand = lambda T: T*np.exp(-T/tau[i])
    integral  = quad(integrand, 0, msa[i])

    A[i] = mass[i] / integral[0]

sfr = A * msa*np.exp(-msa/tau)




# =============================================================================
# PLOT 1
# =============================================================================


cb = tau


xlin = np.array([7, 12])
ylin = slope*xlin + intercept

plt.figure(figsize=(12, 10))
plt.scatter(np.log10(mass), np.log10(sfr), s=5, c=cb) # plot templates
plt.plot(xlin, ylin, color='k') # plot main sequence straight line

plt.xlim(7, 12)
plt.ylim(-1, 3)

plt.colorbar()
plt.show()

# =============================================================================
# PLOT 2
# =============================================================================

plt.figure(figsize=(12, 10))
plt.scatter(np.log10(mass), np.log10(sfr), s=5, c=cb) # plot templates
plt.plot(xlin, ylin, color='k') # plot main sequence straight line

#plt.xlim(7, 12)
#plt.ylim(-1, 3)

plt.show()


# =============================================================================
#np.save('msa', msa)
#np.save('tau', tau)
#np.save('mass', mass)
#np.save('A', A)
#np.save('sfr', sfr)
# =============================================================================

# =============================================================================
# main sequence determined from speagle et al. 2014
# =============================================================================

# original
#intercept = -8.
#slope = 1.

slope2 = 0.84 - 0.026*ageUniv*1e-9 # 0.7560501633562806
intercept2 = -(6.51 - 0.11*ageUniv*1e-9) # -6.1548276141996485

slope3 = 0.84 - 0.026*ageUniv5*1e-9 # 0.7560501633562806
intercept3 = -(6.51 - 0.11*ageUniv5*1e-9) # -6.1548276141996485

xlin = np.array([7, 12])
ylin = slope*xlin + intercept
ylin2 = slope2*xlin + intercept2
ylin3 = slope3*xlin + intercept3

plt.figure(figsize=(12, 10))
plt.plot(xlin, ylin, color='k', label='original') # plot original MS
plt.plot(xlin, ylin2, color='r', label='Speagle, z=2') # plot speagle MS @ z=2
plt.plot(xlin, ylin3, color='b', label='Speagle, z=5') # plot speagle MS @ z=5

plt.xlim(7, 12)
plt.ylim(-1, 3)

plt.legend()
plt.show()

print(ageUniv)
print(ageUniv5)

print(slope2, intercept2)














