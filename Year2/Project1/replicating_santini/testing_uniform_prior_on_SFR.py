#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 09:45:27 2020

@author: lester
"""

# =============================================================================
# DE SFR - uniform prior
# =============================================================================


import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc

nObj = 10

#intercept = -8.
#slope = 1.
#scatter = 0.3

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
cosmo = cd.set_omega_k_0(cosmo)
zObs=2.0
ageUniv = cd.age(zObs, **cosmo)/cc.yr_s
#ageUniv5 = cd.age(5.0, **cosmo)/cc.yr_s

msa = 10**np.random.uniform(low=6, high=np.log10(ageUniv), size=nObj)                      # yrs
#tau = 10**np.random.uniform(low=7, high=10.5, size=nObj)                      # yrs
mass = 10**np.random.uniform(low=7, high=12, size=nObj)                         # solar masses

sfr = 10**np.random.uniform(low=-2, high=2, size=nObj)  # solar masses per year



print(msa[0], mass[0], sfr[0])


tau_test = 10**np.linspace(5, 11, 1000)


# =============================================================================
# plots
# =============================================================================
def calc_sfr(tau, i):
    return ((mass[i]*msa[i]*np.exp(-msa[i]/tau)) / (tau*(tau-(np.exp(-msa[i]/tau)*(tau+msa[i])))))


sfr_min = np.log10(((mass*msa*np.exp(-msa/min(tau_test))) / (min(tau_test)*(min(tau_test)-(np.exp(-msa/min(tau_test))*(min(tau_test)+msa))))))

sfr_max = np.log10(((mass*msa*np.exp(-msa/max(tau_test))) / (max(tau_test)*(max(tau_test)-(np.exp(-msa/max(tau_test))*(max(tau_test)+msa))))))

print(sfr_min)
print(sfr_max)

idx = (sfr_min<0) & (sfr_max>0)
print(idx)

msa = msa[idx]
mass = mass[idx]
sfr = sfr[idx]


for i in range(len(msa)):
    plt.plot(np.log10(tau_test), np.log10(calc_sfr(tau_test, i)))
plt.show()

for i in range(len(msa)):
    plt.plot(np.log10(tau_test), np.log10(calc_sfr(tau_test, i)))
    plt.ylim(-1, 1)
plt.show()


# =============================================================================
# calculate log tau
# =============================================================================

tau = []

for i in range(len(msa)):
    def f(tau):
#        return ((mass[i]*msa[i]*np.exp(-msa[i]/tau)) / (tau*(tau-(np.exp(-msa[i]/tau)*(tau+msa[i]))))) - sfr[i]
        return np.log10(((mass[i]*msa[i]*np.exp(-msa[i]/tau)) / (tau*(tau-(np.exp(-msa[i]/tau)*(tau+msa[i])))))) - np.log10(sfr[i])


    tau.append(fsolve(f, 1e8)[0])

tau = np.array(tau)
#print(tau)

plt.hist(np.log10(tau))
plt.show()

def calc_sfr2(tau):
    return ((mass*msa*np.exp(-msa/tau)) / (tau*(tau-(np.exp(-msa/tau)*(tau+msa)))))


idx2 = abs(np.log10(calc_sfr2(tau)))<2.0

plt.hist(np.log10(tau)[idx2])
plt.show()

plt.hist(np.log10(sfr), alpha=0.3)
plt.hist(np.log10(calc_sfr2(tau))[idx2], alpha=0.3)
plt.show()

plt.scatter(np.log10(sfr)[idx2], np.log10(calc_sfr2(tau))[idx2])
plt.show()

plt.hexbin(np.log10(sfr)[idx2], np.log10(calc_sfr2(tau))[idx2])
plt.show()





# =============================================================================
# different approach, sampling mass, tau, msa then calculating sfr, but varying initial samples
# =============================================================================
#%%

msa = 10**np.random.uniform(low=6, high=np.log10(ageUniv), size=nObj)                      # yrs
tau = 10**np.random.uniform(low=6, high=7, size=nObj)                      # yrs
mass = 10**np.random.uniform(low=7, high=12, size=nObj)                         # solar masses
sfr = ((mass*msa*np.exp(-msa/tau)) / (tau*(tau-(np.exp(-msa/tau)*(tau+msa)))))

plt.scatter(np.log10(mass), np.log10(sfr), alpha=0.1)
plt.xlim(7, 12)
plt.ylim(-3, 3)
plt.show()












