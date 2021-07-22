#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:52:45 2020

@author: lester
"""


import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D




# =============================================================================
# the model (redshift dependent alpha, NO mass dependent scatter (k=1))
# =============================================================================
'''
n = 100
mass = np.random.normal(9, 0.5, n)
z = np.random.uniform(1, 6, n)

alpha_a = -9.0
alpha_b = 0.5

beta_a = 1.0 #0.7
beta_b = 0.0 #0.1

alpha = alpha_a + alpha_b*z
beta = beta_a + beta_b*z
sfr = alpha + mass*beta + np.random.normal(0, 0.3, n)
'''

# =============================================================================
# the model (HOGG + redshift dependent alpha, NO mass dependent scatter (k=1))
# =============================================================================

n = 200
mass = np.random.normal(9, 0.5, n)
z = np.random.uniform(1, 6, n)

alpha_a = -9.0
alpha_b = 0.5

beta_a = 1.0 #0.7
beta_b = 0.0 #0.1

alpha = alpha_a + alpha_b*z
beta = beta_a + beta_b*z

sig0 = 0.3

pbad = 0.3
nBad = np.random.binomial(n,pbad)
outlier_mean = 0
outlier_sigma = 2

tempIdx = np.random.choice(np.fromiter(range(n),np.int),size=nBad,replace=False)
good = np.ones(n,np.bool)
good[tempIdx] = False

sfr = np.zeros_like(mass)
sfr[good] = alpha[good] + mass[good]*beta[good] + np.random.normal(0, sig0, sum(good))
sfr[~good] = np.random.normal(outlier_mean, outlier_sigma, nBad)



# =============================================================================
# reading in scenario 7 for mass distribution
# =============================================================================
'''
options = ['z1p3', 'z2p0', 'z3p0', 'z4p0', 'z5p0']
##options = ['z1p3']
#options = ['z1p3','z5p0']

data = {}
for i in range(len(options)):
    data_temp = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data/scenario_7_7_data_{}.p'.format(options[i]),'r'))
    plt.hist(data_temp['mass_BEAGLE_stellar'], label=options[i], histtype=u'step')

    for key in data_temp.keys():
        
        if i == 0:
            data[key] = []
            data[key].append(data_temp[key])
        else:
            data[key].append(data_temp[key])
            
for key in data.keys():          
    data[key] = np.concatenate(data[key])

plt.legend()
plt.show()    

print(data.keys())
plt.title('Scenario 7 7 TOTAL')
plt.hist(data['mass_BEAGLE_stellar'])
plt.show()

# =============================================================================
# the model with masses and redshifts from scenario 7 7 
# =============================================================================

n = len(data['mass_BEAGLE_stellar'])
mass = data['mass_BEAGLE_stellar']
z = data['redshift_BEAGLE']

alpha_a = -9.0
alpha_b = 0.5

beta_a = 1.0 #0.7
beta_b = 0.0 #0.1

alpha = alpha_a + alpha_b*z
beta = beta_a + beta_b*z
sfr = alpha + mass*beta + np.random.normal(0, 0.3, n)
'''
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
x_GMM_3d = np.array([mass]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
y_GMM_3d = np.array([sfr]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
z_GMM_3d = np.array([z]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))

xsig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.01)) #+ np.random.normal(0, 0.02, (n,3)))
ysig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.01)) #+ np.random.normal(0, 0.02, (n,3)))
zsig_GMM_3d = abs(np.full_like(x_GMM_3d, 0.01)) #+ np.random.normal(0, 0.02, (n,3)))

xycov_GMM_3d = np.full_like(x_GMM_3d, 0.0) #+ np.random.uniform(0, 0, (n,3))
xzcov_GMM_3d = np.full_like(x_GMM_3d, 0.0) #+ np.random.uniform(0, 0, (n,3))
yzcov_GMM_3d = np.full_like(x_GMM_3d, 0.0) #+ np.random.uniform(0, 0, (n,3))

amp_GMM_3d = np.random.uniform(0, 1, (n,3))
amp_GMM_3d_norm = np.sum(amp_GMM_3d, axis=1)
amp_GMM_3d = amp_GMM_3d/amp_GMM_3d_norm[:, np.newaxis]

data = {'x_GMM_3d':x_GMM_3d, 'y_GMM_3d':y_GMM_3d, 'z_GMM_3d':z_GMM_3d, 'xsig_GMM_3d':xsig_GMM_3d, 'ysig_GMM_3d':ysig_GMM_3d, 'zsig_GMM_3d':zsig_GMM_3d, 'xycov_GMM_3d':xycov_GMM_3d, 'xzcov_GMM_3d':xzcov_GMM_3d, 'yzcov_GMM_3d':yzcov_GMM_3d, 'amp_GMM_3d':amp_GMM_3d}

#pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data/mock_hogg_redshift_100.p','w'))

# mock_000
#006 was same without scatter on sfr

#009 was with real mass and redshifts
# 101 - 110 are to test for bias

# mock_hogg_redshift_001 
# 001 is to get the code working...
# 100-109 are for bias checks

test = np.polyfit(mass, sfr, 1)
print(test)




