#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 21:58:26 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

# =============================================================================
# load the templates
# =============================================================================

msa = np.load('msa.npy')
tau = np.load('tau.npy')
tau_exp = np.load('tau_exp.npy')
mass = np.load('mass.npy')
A = np.load('A.npy')
sfr = np.load('sfr.npy')
closest_ind = np.load('closest_ind.npy')

subset = sfr >= 0

msa = msa[subset]
tau = tau[subset]
tau_exp = tau_exp[subset]
mass = mass[subset]
A = A[subset]
sfr = sfr[subset]

msa = msa[closest_ind]
tau = tau[closest_ind]
tau_exp = tau_exp[closest_ind]
mass = mass[closest_ind]
A = A[closest_ind]
sfr = sfr[closest_ind]

# =============================================================================
# plots of how each parameter is sampled
# =============================================================================

plt.title('main sequence')
plt.scatter(np.log10(mass), np.log10(sfr))
plt.show()

plt.title('log msa')
plt.hist(np.log10(msa))
plt.show()

plt.title('log tau')
plt.hist(np.log10(tau))
plt.show()

plt.title('tauexp')
plt.hist(np.log10(tau_exp))
plt.show()


# =============================================================================
# output to file
# =============================================================================

outputDict = {}

outputDict['id']                = np.array(range(len(mass)))+1
outputDict['max_stellar_age']   = np.log10(msa)
outputDict['tau']               = np.log10(tau)
outputDict['tau_exp']           = np.log10(tau_exp)
outputDict['metallicity']       = np.random.uniform(-2.1, 0.3, len(mass))
outputDict['mass']              = np.log10(mass)

outputDict['nebular_logU']      = np.random.uniform(-4, -1, len(mass))
outputDict['nebular_xi']        = np.random.uniform(0.1, 0.5, len(mass))

outputDict['tauV_eff']          = np.random.uniform(0, 2, len(mass))

outputDict['A']                 = A
outputDict['sfr']               = sfr
outputDict['closest_ind']       = closest_ind

outputTable = Table(outputDict)
#outputTable.write("mock_MS_parameters_100_LE.fits", overwrite=True)



























