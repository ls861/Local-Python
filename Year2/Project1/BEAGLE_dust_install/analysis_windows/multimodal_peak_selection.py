#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:28:49 2021

@author: lester
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm,truncnorm,multivariate_normal
import time
from astropy.io import fits


x = np.linspace(-25, 25, 1000)

a1 = 5
b1 = 3
c1 = 1

a2 = 10
b2 = 4
c2 = 2

y1 = a1 * np.exp(-((x-b1)**2)/(2*c1**2))
y2 = a2 * np.exp(-((x-b2)**2)/(2*c2**2))


plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y1*y2)
plt.plot(x, (a1*c1*np.sqrt(2*np.pi))*norm.pdf(x, b1, c1))
plt.show()

#plt.plot(x, norm.pdf(x, b1, c1))
#plt.plot(x, norm.pdf(x, b2, c2))
#plt.plot(x, norm.pdf(x, b1, c1) * norm.pdf(x, b2, c2))
#plt.show()




#%%


# =============================================================================
# coordinate grid
# =============================================================================

x = np.linspace(0, 20, 20)
y = np.linspace(-10, 10, 20)
z = np.linspace(0, 10, 20)

g = np.meshgrid(x,y,z, sparse=False, indexing='ij')

coords = np.empty([0, 3])

for i in range(len(x)):
    for j in range(len(y)):
        for k in range(len(z)):
            coord = np.array([[g[0][i][j][k], g[1][i][j][k], g[2][i][j][k]]])
            coords = np.concatenate((coords, coord))

# =============================================================================
# gmm pdfs
# =============================================================================

n_gmm = 1.0
pi_gmm = 1.0
mean = np.array([10, 0, 5])
cov = np.matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
pdf_gmm = multivariate_normal.pdf(coords, mean, cov)

# =============================================================================
# ms pdfs
# =============================================================================

n_ms = 1.0
alpha = 1.0
beta = 0.8
sigma = 0.3
pdf_ms = norm.pdf(coords[:,1], beta*(coords[:,0]-9.7) + alpha, sigma)

# =============================================================================
# mass pdfs
# =============================================================================

n_m = 1.0
mass_mu = np.array([9.0, 9.0, 9.0])
mass_sigma = np.array([0.5, 0.5, 0.5])
mass_pi = np.array([0.3, 0.3, 0.4])

pdf_m = mass_pi[0] * norm.pdf(coords[:,0], mass_mu[0], mass_sigma[0]) + \
        mass_pi[1] * norm.pdf(coords[:,0], mass_mu[1], mass_sigma[1]) + \
        mass_pi[2] * norm.pdf(coords[:,0], mass_mu[2], mass_sigma[2])
    
# =============================================================================
# total
# =============================================================================

pdf = (n_gmm * n_ms * n_m * pi_gmm * pdf_gmm * pdf_ms * pdf_m) # volume of cube **3 not needed, can be factored into nomalisations if necessary

print(sum(pdf))


t0 = time.time()
t1 = time.time()



#%%
# =============================================================================
# real stuff
# =============================================================================

def pdf(s, idx, G, n):
    
    x = np.linspace(0, 20, n)
    y = np.linspace(-10, 10, n)
    z = np.linspace(0, 10, n)
    g = np.meshgrid(x,y,z, sparse=False, indexing='ij')
    coords = np.empty([0, 3])
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                coord = np.array([[g[0][i][j][k], g[1][i][j][k], g[2][i][j][k]]])
                coords = np.concatenate((coords, coord))
    print(coords)
    # =============================================================================
    # gmm pdfs
    # =============================================================================
    n_gmm = 1.0
    pi_gmm = s['amp_GMM_3d'][idx,G]
    mean = np.array([s['x_GMM_3d'][idx,G],s['y_GMM_3d'][idx,G],s['z_GMM_3d'][idx,G]])
    cov = np.array([[np.power(s['xsig_GMM_3d'][idx,G],2), s['xycov_GMM_3d'][idx,G], s['xzcov_GMM_3d'][idx,G]],[s['xycov_GMM_3d'][idx,G], np.power(s['ysig_GMM_3d'][idx,G],2), s['yzcov_GMM_3d'][idx,G]],[s['xzcov_GMM_3d'][idx,G], s['yzcov_GMM_3d'][idx,G], np.power(s['zsig_GMM_3d'][idx,G],2)]])
    pdf_gmm = multivariate_normal.pdf(coords, mean, cov)
#    print('pdf_gmm', sum(pdf_gmm))
    
    # =============================================================================
    # ms pdfs
    # =============================================================================
    n_ms = 1.0
    alpha = np.full(np.shape(coords[:,0]), 0.0)
    beta = np.full(np.shape(coords[:,0]), 0.0)
    sigma = np.full(np.shape(coords[:,0]), 0.0)

    # ms per z bin
    alpha = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.98, alpha)
    beta = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.86, beta)
    sigma = np.where((coords[:,2]>-1.0)&(coords[:,2]<2.0), 0.31, sigma)

    alpha = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 1.02, alpha)
    beta = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 0.89, beta)
    sigma = np.where((coords[:,2]>2.0)&(coords[:,2]<3.0), 0.21, sigma)

    alpha = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 1.10, alpha)
    beta = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 0.72, beta)
    sigma = np.where((coords[:,2]>3.0)&(coords[:,2]<4.0), 0.22, sigma)

    alpha = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 1.68, alpha)
    beta = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 0.93, beta)
    sigma = np.where((coords[:,2]>3.0)&(coords[:,2]<5.0), 0.33, sigma)

    alpha = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 2.67, alpha)
    beta = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 1.54, beta)
    sigma = np.where((coords[:,2]>4.0)&(coords[:,2]<99.0), 0.46, sigma)
        
    pdf_ms = norm.pdf(coords[:,1], beta*(coords[:,0]-9.7) + alpha, sigma)
#    print('pdf_ms', sum(pdf_ms))
    
    # =============================================================================
    # mass pdfs
    # =============================================================================
    n_m = 1.0
    mass_mu = np.full((len(alpha),3), 0.0)
    mass_sigma = np.full((len(alpha),3), 0.0)
    mass_pi = np.full((len(alpha),3), 0.0)

    # mass per z bin
    mass_mu[:,0] = np.where(coords[:,2]<99, 8.0, mass_mu[:,0])
    mass_mu[:,1] = np.where(coords[:,2]<99, 8.0, mass_mu[:,1])
    mass_mu[:,2] = np.where(coords[:,2]<99, 8.0, mass_mu[:,2])
    mass_sigma[:,0] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,0])
    mass_sigma[:,1] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,1])
    mass_sigma[:,2] = np.where(coords[:,2]<99, 0.5, mass_sigma[:,2])    
    mass_pi[:,0] = np.where(coords[:,2]<99, 0.3, mass_pi[:,0])
    mass_pi[:,1] = np.where(coords[:,2]<99, 0.3, mass_pi[:,1])
    mass_pi[:,2] = np.where(coords[:,2]<99, 0.3, mass_pi[:,2])

    pdf_m = mass_pi[:,0] * norm.pdf(coords[:,0], mass_mu[:,0], mass_sigma[:,0]) + \
            mass_pi[:,1] * norm.pdf(coords[:,0], mass_mu[:,1], mass_sigma[:,1]) + \
            mass_pi[:,2] * norm.pdf(coords[:,0], mass_mu[:,2], mass_sigma[:,2])
            
    pdf_m = 1.0

    # =============================================================================
    # Total
    # =============================================================================
    pdf = (n_gmm * n_ms * n_m * pi_gmm * pdf_gmm * pdf_ms * pdf_m) # volume of cube **3 not needed, can be factored into nomalisations if necessary
    return sum(pdf)



fileName = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_29_clusters_z1p25-2p0.fits'
s = fits.open(fileName)[1].data



distances = []

for i in range(len(s)):
#for i in [264]:    
#    print(s['field_AD'][i], s['id_AD'][i], i)    
#    
#    d1 = np.sqrt(    (s['x_GMM_3d'][i][0]-s['x_GMM_3d'][i][1])**2 + \
#                     (s['y_GMM_3d'][i][0]-s['y_GMM_3d'][i][1])**2 + \
#                     (s['z_GMM_3d'][i][0]-s['z_GMM_3d'][i][1])**2 )
#
#    d2 = np.sqrt(    (s['x_GMM_3d'][i][0]-s['x_GMM_3d'][i][2])**2 + \
#                     (s['y_GMM_3d'][i][0]-s['y_GMM_3d'][i][2])**2 + \
#                     (s['z_GMM_3d'][i][0]-s['z_GMM_3d'][i][2])**2 )
#
#    d3 = np.sqrt(    (s['x_GMM_3d'][i][1]-s['x_GMM_3d'][i][2])**2 + \
#                     (s['y_GMM_3d'][i][1]-s['y_GMM_3d'][i][2])**2 + \
#                     (s['z_GMM_3d'][i][1]-s['z_GMM_3d'][i][2])**2 )    
#  
#    distances.append(max(d1, d2, d3))
    
    for G in range(3):
        print(pdf(s, i, G, 30))
    print(' ')

#plt.hist(distances, bins=50)
#plt.show()





#%%
# =============================================================================
# kelly code normalisation issue
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm,truncnorm,multivariate_normal
import time
from astropy.io import fits

# plot normalised gaussians
x = np.linspace(-25, 25, 1000)
x_scatter = np.linspace(0, 10, 10)

mu1 = 3.0
sig1 = 4.0
n1 = 1.0 / (sig1 * np.sqrt(2*np.pi))

mu2 = 10.0
sig2 = 7.0
n2 = 1.0 / (sig2 * np.sqrt(2*np.pi))

y1 = n1 * np.exp(-((x-mu1)**2)/(2*sig1**2))
y2 = n2 * np.exp(-((x-mu2)**2)/(2*sig2**2))

mu12 = (mu1*(sig2**2) + mu2*(sig1**2)) / ((sig1**2)+(sig2**2))
sig12 = np.sqrt(((sig1**2)*(sig2**2))  / ((sig1**2)+(sig2**2)))
n12 = 1.0 / (sig12 * np.sqrt(2*np.pi))
y12 = n12 * np.exp(-((x-mu12)**2)/(2*sig12**2))


S_12 = (1.0/np.sqrt(2*np.pi*((sig1**2)+(sig2**2)))) * np.exp(-((mu1-mu2)**2)/(2*((sig1**2)+(sig2**2))))

N = (2*np.pi*sig1*sig2)/((np.sqrt(2*np.pi))*np.sqrt(((sig1**2)*(sig2**2))  / ((sig1**2)+(sig2**2))))
N = 1.0 / (n1*sig1*np.sqrt(2*np.pi))


plt.figure(figsize=(10,10))
plt.plot(x, y1, label='pdf1')
plt.plot(x, y2, label='pdf2')
plt.plot(x, y12, label='pdf12')
plt.plot(x, y1*y2/S_12, label='1x2 / S')
plt.plot(x, y1*y2, label='1x2')

plt.scatter(x_scatter, norm.pdf(x_scatter, mu1, sig1))
plt.scatter(x_scatter, norm.pdf(x_scatter, mu2, sig2))
plt.scatter(x_scatter, norm.pdf(x_scatter, mu12, sig12))
#plt.scatter(x_scatter, norm.pdf(x_scatter, mu1, sig1)*norm.pdf(x_scatter, mu2, sig2), marker='x')
#plt.scatter(x_scatter, N*norm.pdf(x_scatter, mu1, sig1)*norm.pdf(x_scatter, mu2, sig2), marker='x')

plt.legend()
plt.show()


# =============================================================================
# 3d
# =============================================================================

mu3 = -5.0
sig3 = 14.0
n3 = 1.0 / (sig3 * np.sqrt(2*np.pi))
y3 = n3 * np.exp(-((x-mu3)**2)/(2*sig3**2))

mu123 = (mu12*(sig3**2) + mu3*(sig12**2)) / ((sig12**2)+(sig3**2))
sig123 = np.sqrt(((sig12**2)*(sig3**2))  / ((sig12**2)+(sig3**2)))
n123 = 1.0 / (sig123 * np.sqrt(2*np.pi))
y123 = n123 * np.exp(-((x-mu123)**2)/(2*sig123**2))

S_123 = (1.0/np.sqrt(2*np.pi*((sig12**2)+(sig3**2)))) * np.exp(-((mu12-mu3)**2)/(2*((sig12**2)+(sig3**2))))


plt.figure(figsize=(10,10))
plt.plot(x, y1, label='pdf1')
plt.plot(x, y2, label='pdf2')
plt.plot(x, y3, label='pdf3')
plt.plot(x, y123, label='pdf123')

plt.plot(x, y1*y2*y3/(S_12*S_123), label='1x2x3', linestyle='dashed')
plt.plot(x, y1*y2*y3, label='1x2x3 uncorrected', linestyle='dashed')

plt.scatter(x_scatter, norm.pdf(x_scatter, mu1, sig1))
plt.scatter(x_scatter, norm.pdf(x_scatter, mu2, sig2))
plt.scatter(x_scatter, norm.pdf(x_scatter, mu3, sig3))
plt.scatter(x_scatter, norm.pdf(x_scatter, mu123, sig123))
#plt.scatter(x_scatter, norm.pdf(x_scatter, mu1, sig1)*norm.pdf(x_scatter, mu2, sig2), marker='x')
#plt.scatter(x_scatter, N*norm.pdf(x_scatter, mu1, sig1)*norm.pdf(x_scatter, mu2, sig2), marker='x')

plt.legend()
plt.show()














print(np.log(norm.pdf(3, 4, 2)))
print(norm.logpdf(3, 4, 2))



print(np.exp(np.log(norm.pdf(3, 4, 2))))
print(norm.pdf(3, 4, 2))











