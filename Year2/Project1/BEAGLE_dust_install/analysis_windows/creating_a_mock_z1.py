#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: lester
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


for m in np.arange(10):
    
    # =============================================================================
    # the model (HOGG + redshift dependent alpha, NO mass dependent scatter (k=1))
    # values taken from full scenario 29 fit (ssfr alpha z1.25 - 6.0)
    # =============================================================================
    
    n = 100
    mass = np.random.normal(8.5, 0.5, n)
    z = np.random.uniform(1.25, 2.0, n)
    
    alpha_a = 0.9
    alpha_b = 0.0
    
    beta_a = 1.0
    beta_b = 0.0
    
    alpha = alpha_a + alpha_b*z
    beta = beta_a + beta_b*z
    
    sig0 = 0.25
    
    pbad = 0.2
    nBad = np.random.binomial(n,pbad)
    outlier_mean = -2.0
    outlier_sigma = 1.2
    
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
    # kelly input
    # =============================================================================
    x_GMM_3d = np.array([mass]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
    y_GMM_3d = np.array([sfr]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
    z_GMM_3d = np.array([z]*3).transpose() #+ np.random.normal(0, 0.1, (n,3))
    
    xsig_GMM_3d = np.random.normal(0.15, 0.05, (n,3))
    xsig_GMM_3d = np.where(xsig_GMM_3d<0.01, 0.01, xsig_GMM_3d)
    ysig_GMM_3d = np.random.normal(0.2, 0.05, (n,3))
    ysig_GMM_3d = np.where(ysig_GMM_3d<0.01, 0.01, ysig_GMM_3d)
    zsig_GMM_3d = np.random.normal(0.07, 0.05, (n,3))
    zsig_GMM_3d = np.where(zsig_GMM_3d<0.01, 0.01, zsig_GMM_3d)
    
    xycov_GMM_3d = np.random.normal(-0.3, 0.4, (n,3))
    xycov_GMM_3d = np.where(np.abs(xycov_GMM_3d)>1.0, -0.3, xycov_GMM_3d) * xsig_GMM_3d * ysig_GMM_3d
    xzcov_GMM_3d = np.random.normal(0.3, 0.4, (n,3))
    xzcov_GMM_3d = np.where(np.abs(xzcov_GMM_3d)>1.0, 0.3, xzcov_GMM_3d) * xsig_GMM_3d * zsig_GMM_3d
    yzcov_GMM_3d = np.random.normal(0.3, 0.4, (n,3))
    yzcov_GMM_3d = np.where(np.abs(yzcov_GMM_3d)>1.0, 0.3, yzcov_GMM_3d) * ysig_GMM_3d * zsig_GMM_3d
    
    amp_GMM_3d = np.random.uniform(0, 1, (n,3))
    amp_GMM_3d_norm = np.sum(amp_GMM_3d, axis=1)
    amp_GMM_3d = amp_GMM_3d/amp_GMM_3d_norm[:, np.newaxis]
    
    
    # =============================================================================
    # sorting out the eigenvalue issue
    # =============================================================================
    count = 0
    for i in range(len(x_GMM_3d)):
        for j in np.arange(3):
            
            cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                      [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                      [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]
                
            while not np.all(np.linalg.eigvals(cov) > 0):
                count+=1
                print(count)
                
                xycov_GMM_3d[i][j] = np.random.normal(-0.3, 0.4)
                xycov_GMM_3d[i][j] = np.where(np.abs(xycov_GMM_3d[i][j])>1.0, -0.3, xycov_GMM_3d[i][j]) * xsig_GMM_3d[i][j] * ysig_GMM_3d[i][j]
                xzcov_GMM_3d[i][j] = np.random.normal(0.3, 0.4)
                xzcov_GMM_3d[i][j] = np.where(np.abs(xzcov_GMM_3d[i][j])>1.0, 0.3, xzcov_GMM_3d[i][j]) * xsig_GMM_3d[i][j] * zsig_GMM_3d[i][j]
                yzcov_GMM_3d[i][j] = np.random.normal(0.3, 0.4)
                yzcov_GMM_3d[i][j] = np.where(np.abs(yzcov_GMM_3d[i][j])>1.0, 0.3, yzcov_GMM_3d[i][j]) * ysig_GMM_3d[i][j] * zsig_GMM_3d[i][j]                
                
                cov = [[xsig_GMM_3d[i][j]**2,xycov_GMM_3d[i][j],xzcov_GMM_3d[i][j]],\
                      [xycov_GMM_3d[i][j],ysig_GMM_3d[i][j]**2,yzcov_GMM_3d[i][j]],\
                      [xzcov_GMM_3d[i][j],yzcov_GMM_3d[i][j],zsig_GMM_3d[i][j]**2]]
    
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
    
    # =============================================================================
    # adding a double peaked object
    # =============================================================================
    #
    #mass1 = np.random.normal(9, 0.5)
    #z1 = 1.625
    #alpha1 = np.log10(alpha_a*(1+z1)**alpha_b) - 9.0 + 9.7
    #beta1 = beta_a + beta_b*z1
    #sfr1 = alpha1 + (mass1-9.7)*beta1 + np.random.normal(0, sig0)
    #
    #mass2 = np.random.normal(9, 0.5)
    #z2 = 4.5
    #alpha2 = np.log10(alpha_a*(1+z2)**alpha_b) - 9.0 + 9.7
    #beta2 = beta_a + beta_b*z2
    #sfr2 = alpha2 + (mass2-9.7)*beta2 + np.random.normal(0, sig0)
    #
    #x_GMM_3d[0] = np.array([mass1, mass1, mass2]).transpose() + np.random.normal(0, 0.1, 3)
    #y_GMM_3d[0] = np.array([sfr1, sfr1, sfr2]).transpose() + np.random.normal(0, 0.1, 3)
    #z_GMM_3d[0] = np.array([z1, z1, z2]).transpose() + np.random.normal(0, 0.1, 3)
    #amp_GMM_3d[0] = np.array([0.3, 0.3, 0.4])
    
    data = {'x_GMM_3d':x_GMM_3d, 'y_GMM_3d':y_GMM_3d, 'z_GMM_3d':z_GMM_3d, 'xsig_GMM_3d':xsig_GMM_3d, 'ysig_GMM_3d':ysig_GMM_3d, 'zsig_GMM_3d':zsig_GMM_3d, 'xycov_GMM_3d':xycov_GMM_3d, 'xzcov_GMM_3d':xzcov_GMM_3d, 'yzcov_GMM_3d':yzcov_GMM_3d, 'amp_GMM_3d':amp_GMM_3d}
    
    #pickle.dump(data, open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/mock_z1_001.p','w'))
    
    #WINDOWS
    # pickle.dump(data, open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/mock_z1_0{}.p'.format(str(m+11)),'wb'))
    
    
    print(x_GMM_3d[0])
    print(y_GMM_3d[0])
    print(z_GMM_3d[0])
    print(amp_GMM_3d[0])
    print(nBad)
    
    #test = np.polyfit(mass, sfr, 1)
    #print(test)




