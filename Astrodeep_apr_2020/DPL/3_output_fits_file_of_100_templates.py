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

alpha = np.load('alpha.npy')
beta = np.load('beta.npy')
tau = np.load('tau.npy')
mass = np.load('mass.npy')
A = np.load('A.npy')
sfr = np.load('sfr.npy')
closest_ind = np.load('closest_ind.npy')

subset = sfr >= 0

alpha = alpha[subset]
beta = beta[subset]
tau = tau[subset]
mass = mass[subset]
A = A[subset]
sfr = sfr[subset]

alpha = alpha[closest_ind]
beta = beta[closest_ind]
tau = tau[closest_ind]
mass = mass[closest_ind]
A = A[closest_ind]
sfr = sfr[closest_ind]

# =============================================================================
# plots of how each parameter is sampled
# =============================================================================

plt.title('main sequence')
plt.scatter(np.log10(mass), np.log10(sfr))
plt.show()

plt.title('log alpha')
plt.hist(np.log10(alpha))
plt.show()

plt.title('log beta')
plt.hist(np.log10(beta))
plt.show()

plt.title('tau')
plt.hist(tau)
plt.show()


# =============================================================================
# output to file
# =============================================================================

outputDict = {}

outputDict['id']                = np.array(range(len(mass)))+1
outputDict['tau']               = np.log10(tau)
outputDict['dpl_alpha']         = alpha
outputDict['dpl_beta']          = beta
outputDict['metallicity']       = np.random.uniform(-2.1, 0.3, len(mass))
outputDict['mass']              = np.log10(mass)

outputDict['nebular_logU']      = np.random.uniform(-4, -1, len(mass))
outputDict['nebular_xi']        = np.random.uniform(0.1, 0.5, len(mass))

outputDict['tauV_eff']          = np.random.uniform(0, 2, len(mass))

outputDict['A']                 = A
outputDict['sfr']               = sfr
outputDict['closest_ind']       = closest_ind

outputTable = Table(outputDict)
#outputTable.write("mock_MS_parameters_100_DPL.fits", overwrite=True)



print(sfr)























