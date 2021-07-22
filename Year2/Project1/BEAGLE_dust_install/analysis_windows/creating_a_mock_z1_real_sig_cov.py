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

m=0

n = 400
mass = np.random.normal(8.5, 0.7, n)
z = np.random.uniform(1.25, 2.0, n)

alpha_a = 1.0
alpha_b = 0.0

beta_a = 1.0
beta_b = 0.0

alpha = alpha_a + alpha_b*z
beta = beta_a + beta_b*z

sig0 = 0.25

pbad = 0.2
nBad = np.random.binomial(n,pbad)
outlier_mean = -3.0
outlier_sigma = 2.0

tempIdx = np.random.choice(np.fromiter(range(n),np.int),size=nBad,replace=False)
good = np.ones(n,np.bool)
good[tempIdx] = False

sfr = np.zeros_like(mass)
sfr[good] = alpha[good] + (mass[good]-9.7)*beta[good] + np.random.normal(0, sig0, sum(good))
sfr[~good] = np.random.normal(outlier_mean, outlier_sigma, nBad)

plt.title(m+11)
plt.hist(mass)
plt.show()
plt.title(m+11)
plt.hist(sfr)
plt.show()
plt.title(m+11)
plt.hist(z)
plt.show()

# =============================================================================
# plot
# =============================================================================
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(mass, z, sfr, c=z)
plt.show()

plt.title(m+11)
plt.scatter(mass, sfr, c=z)
plt.colorbar()
plt.show()


# =============================================================================
# find closest real object from z1 bin
# =============================================================================

#WINDOWS
scenarioA = '29'
field = 'clusters'
z_bin = 'z1p25-2p0'
z_lower = 1.25
z_upper = 2.0
with open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin), 'rb') as f:
    real_data = pickle.load(f, encoding='latin1')

x_real = np.empty(len(real_data['amp_GMM_3d']))
y_real = np.empty(len(real_data['amp_GMM_3d']))
z_real = np.empty(len(real_data['amp_GMM_3d']))

xsig_real = np.empty(len(real_data['amp_GMM_3d']))
ysig_real = np.empty(len(real_data['amp_GMM_3d']))
zsig_real = np.empty(len(real_data['amp_GMM_3d']))
xycov_real = np.empty(len(real_data['amp_GMM_3d']))
xzcov_real = np.empty(len(real_data['amp_GMM_3d']))
yzcov_real = np.empty(len(real_data['amp_GMM_3d']))

# choose the highest probabiliy cloud per real object
for i in range(len(real_data['amp_GMM_3d'])):
    maxind = np.argmax(real_data['amp_GMM_3d'][i])

    x_real[i] = real_data['x_GMM_3d'][i, maxind]
    y_real[i] = real_data['y_GMM_3d'][i, maxind]
    z_real[i] = real_data['z_GMM_3d'][i, maxind]
    
    xsig_real[i] = real_data['xsig_GMM_3d'][i, maxind]
    ysig_real[i] = real_data['ysig_GMM_3d'][i, maxind]
    zsig_real[i] = real_data['zsig_GMM_3d'][i, maxind]
    xycov_real[i] = real_data['xycov_GMM_3d'][i, maxind]
    xzcov_real[i] = real_data['xzcov_GMM_3d'][i, maxind]
    yzcov_real[i] = real_data['yzcov_GMM_3d'][i, maxind]
    
x_real_closest_to_mock = np.empty(n)
y_real_closest_to_mock = np.empty(n)
z_real_closest_to_mock = np.empty(n)

xsig_GMM_3d = np.empty(n)
ysig_GMM_3d = np.empty(n)
zsig_GMM_3d = np.empty(n)
xycov_GMM_3d = np.empty(n)
xzcov_GMM_3d = np.empty(n)
yzcov_GMM_3d = np.empty(n)

# per mock object, decide which real object is closest in 3d space
for i in range(n):
    
    dist = np.sqrt((mass[i]-x_real)**2 + (sfr[i]-y_real)**2 + (z[i]-z_real)**2)    
    x_real_closest_to_mock[i] = x_real[np.argmin(dist)]
    y_real_closest_to_mock[i] = y_real[np.argmin(dist)]
    z_real_closest_to_mock[i] = z_real[np.argmin(dist)]
    
    xsig_GMM_3d[i] = xsig_real[np.argmin(dist)]
    ysig_GMM_3d[i] = ysig_real[np.argmin(dist)]
    zsig_GMM_3d[i] = zsig_real[np.argmin(dist)]
    xycov_GMM_3d[i] = xycov_real[np.argmin(dist)]
    xzcov_GMM_3d[i] = xzcov_real[np.argmin(dist)]
    yzcov_GMM_3d[i] = yzcov_real[np.argmin(dist)]

# plot MS showing mock points and means of chosen closest real clouds
plt.figure(figsize=(10,10))
for i in range(n):
    plt.plot((mass[i], x_real_closest_to_mock[i]),(sfr[i], y_real_closest_to_mock[i]))

plt.scatter(x_real, y_real, color='k', s=1)
plt.show()


#%%
# =============================================================================
# kelly input
# =============================================================================
x_GMM_3d = np.array([mass]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
y_GMM_3d = np.array([sfr]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
z_GMM_3d = np.array([z]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))

xsig_GMM_3d = np.array([xsig_GMM_3d]*3).transpose()
ysig_GMM_3d = np.array([ysig_GMM_3d]*3).transpose()
zsig_GMM_3d = np.array([zsig_GMM_3d]*3).transpose()

xycov_GMM_3d = np.array([xycov_GMM_3d]*3).transpose()
xzcov_GMM_3d = np.array([xzcov_GMM_3d]*3).transpose()
yzcov_GMM_3d = np.array([yzcov_GMM_3d]*3).transpose()

amp_GMM_3d = np.random.uniform(0, 1, (n,3))
amp_GMM_3d_norm = np.sum(amp_GMM_3d, axis=1)
amp_GMM_3d = amp_GMM_3d/amp_GMM_3d_norm[:, np.newaxis]


# =============================================================================
# perturbing means based on covariance matrices
# =============================================================================
for i in range(len(x_GMM_3d)):
    for j in np.arange(3):
        
        cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                  [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                  [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]
        mean = [x_GMM_3d[i,j],y_GMM_3d[i,j],z_GMM_3d[i,j]]
        
        xyz = np.random.multivariate_normal(mean, cov)
        
        x_GMM_3d[i][j] = xyz[0]
        y_GMM_3d[i][j] = xyz[1]
        z_GMM_3d[i][j] = xyz[2]

plt.title(m+11)
plt.scatter(x_GMM_3d.flatten(), y_GMM_3d.flatten(), c=z_GMM_3d.flatten())
plt.colorbar()
plt.show()


data = {'x_GMM_3d':x_GMM_3d, 'y_GMM_3d':y_GMM_3d, 'z_GMM_3d':z_GMM_3d, 'xsig_GMM_3d':xsig_GMM_3d, 'ysig_GMM_3d':ysig_GMM_3d, 'zsig_GMM_3d':zsig_GMM_3d, 'xycov_GMM_3d':xycov_GMM_3d, 'xzcov_GMM_3d':xzcov_GMM_3d, 'yzcov_GMM_3d':yzcov_GMM_3d, 'amp_GMM_3d':amp_GMM_3d}

# pickle.dump(data, open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/mock_z1_real_sig_cov_002.p','wb'))


print(x_GMM_3d[0])
print(y_GMM_3d[0])
print(z_GMM_3d[0])
print(amp_GMM_3d[0])
print(nBad)

#test = np.polyfit(mass, sfr, 1)
#print(test)




