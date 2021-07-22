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

#for i in range(len(s)):
for i in [264]:
    for G in range(3):
        print(pdf(s, i, G, 10))
    print(' ')
