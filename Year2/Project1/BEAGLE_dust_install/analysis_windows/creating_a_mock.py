#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: lester
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# the model (HOGG + redshift dependent alpha, NO mass dependent scatter (k=1))
# values taken from full scenario 29 fit (ssfr alpha z1.25 - 6.0)
# =============================================================================

n = 500
mass = np.random.normal(9, 0.5, n)
z = np.random.uniform(1.25, 6, n)

alpha_a = 0.13
alpha_b = 2.18

beta_a = 0.768
beta_b = 0.0

alpha = np.log10(alpha_a*(1+z)**alpha_b) - 9.0 + 9.7
beta = beta_a + beta_b*z

sig0 = 0.25

pbad = 0.242
nBad = np.random.binomial(n,pbad)
outlier_mean = 1.099
outlier_sigma = 0.857

tempIdx = np.random.choice(np.fromiter(range(n),np.int),size=nBad,replace=False)
good = np.ones(n,np.bool)
good[tempIdx] = False

sfr = np.zeros_like(mass)
sfr[good] = alpha[good] + (mass[good]-9.7)*beta[good] + np.random.normal(0, sig0, sum(good))
sfr[~good] = np.random.normal(outlier_mean, outlier_sigma, nBad)

plt.hist(mass)
plt.show()
plt.hist(sfr)
plt.show()
plt.hist(z)
plt.show()

# =============================================================================
# plot
# =============================================================================
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(mass, z, sfr, c=z)
plt.show()

# =============================================================================
# kelly input
# =============================================================================
x_GMM_3d = np.array([mass]*3).transpose() + np.random.normal(0, 0.1, (n,3))
y_GMM_3d = np.array([sfr]*3).transpose() + np.random.normal(0, 0.1, (n,3))
z_GMM_3d = np.array([z]*3).transpose() + np.random.normal(0, 0.1, (n,3))

xsig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.1))
ysig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.1))
zsig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.1))

xycov_GMM_3d = np.full_like(x_GMM_3d, 0.0)
xzcov_GMM_3d = np.full_like(x_GMM_3d, 0.0)
yzcov_GMM_3d = np.full_like(x_GMM_3d, 0.0)

amp_GMM_3d = np.random.uniform(0, 1, (n,3))
amp_GMM_3d_norm = np.sum(amp_GMM_3d, axis=1)
amp_GMM_3d = amp_GMM_3d/amp_GMM_3d_norm[:, np.newaxis]

# =============================================================================
# adding a double peaked object
# =============================================================================

mass1 = np.random.normal(9, 0.5)
z1 = 1.625
alpha1 = np.log10(alpha_a*(1+z1)**alpha_b) - 9.0 + 9.7
beta1 = beta_a + beta_b*z1
sfr1 = alpha1 + (mass1-9.7)*beta1 + np.random.normal(0, sig0)

mass2 = np.random.normal(9, 0.5)
z2 = 4.5
alpha2 = np.log10(alpha_a*(1+z2)**alpha_b) - 9.0 + 9.7
beta2 = beta_a + beta_b*z2
sfr2 = alpha2 + (mass2-9.7)*beta2 + np.random.normal(0, sig0)

x_GMM_3d[0] = np.array([mass1, mass1, mass2]).transpose() + np.random.normal(0, 0.1, 3)
y_GMM_3d[0] = np.array([sfr1, sfr1, sfr2]).transpose() + np.random.normal(0, 0.1, 3)
z_GMM_3d[0] = np.array([z1, z1, z2]).transpose() + np.random.normal(0, 0.1, 3)
amp_GMM_3d[0] = np.array([0.3, 0.3, 0.4])

data = {'x_GMM_3d':x_GMM_3d, 'y_GMM_3d':y_GMM_3d, 'z_GMM_3d':z_GMM_3d, 'xsig_GMM_3d':xsig_GMM_3d, 'ysig_GMM_3d':ysig_GMM_3d, 'zsig_GMM_3d':zsig_GMM_3d, 'xycov_GMM_3d':xycov_GMM_3d, 'xzcov_GMM_3d':xzcov_GMM_3d, 'yzcov_GMM_3d':yzcov_GMM_3d, 'amp_GMM_3d':amp_GMM_3d}

#pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/mock_003.p','w'))

print(x_GMM_3d[0])
print(y_GMM_3d[0])
print(z_GMM_3d[0])
print(amp_GMM_3d[0])

#test = np.polyfit(mass, sfr, 1)
#print(test)



