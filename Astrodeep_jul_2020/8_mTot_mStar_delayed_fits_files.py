#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:17:54 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

massType = '_mStar'
#massType = '_mTot'

sfrType = '_delayed'
#sfrType = ''

massLim = '8p4'
redshift_center     = 1.65 # 2.5, 3.5, 4.5
redshift_width      = 0.35 # either side
chi2_max            = 9999 # 2.5, 9.5
dimension           = '2d' # 2d or 3d - 3d includes redshift in the GMM fit
comment             = '001'

zLim = str(redshift_center).replace(".","p")
wLim = str(redshift_width).replace(".","p")
chi2Lim = str(chi2_max).replace(".","p")



sbf = './from_cluster/linmix_AD_combined/linmix_npy_files_z{}_w{}_chi{}_dim{}{}{}_{}/'.format(zLim, wLim, chi2Lim, dimension, massType, sfrType, comment)

sbf = '{}{}_'.format(sbf, massLim)

GMMid       = np.load(sbf+'id.npy')  
GMMfield    = np.load(sbf+'field.npy')  
pi_err      = np.load(sbf+'pi_err.npy')         # 3x probability of each posterior gaussian
GMMx        = np.load(sbf+'GMMx.npy')           # 3x posterior means per mass
GMMy        = np.load(sbf+'GMMy.npy')           # 3x posterior means per sfr
GMMxsig     = np.load(sbf+'GMMxsig.npy')        # 3x posterior sigmas per mass
GMMysig     = np.load(sbf+'GMMysig.npy')        # 3x posterior sigmas per sfr
GMMxycov    = np.load(sbf+'GMMxycov.npy')       # 3x posterior covar per mass-sfr pair



print('HEATPLOTS')



samples = 1000
xlow = 8.0
xhigh = 10.5
ylow = -1.5
yhigh = 4
# objs to plot
#objs = np.random.choice(range(len(GMMx)), 2) # random selection
objs = np.array(range(len(GMMx))) # all





idx = np.isin(GMMid, [ 2,  292,   17,   59,  785, 1766, 1800, 1307, 2180, 2372, 2426, 2439, 2488, 2499])
idx = np.isin(GMMid, [292, 17])

print(idx)



for obj in objs[idx]:
#for obj in [0]:
#    print(obj)
    print('ID: ', GMMfield[obj], GMMid[obj])
   
    hp_x = []
    hp_y = []

    draws = np.random.choice([0, 1, 2], samples, p=pi_err[obj]/sum(pi_err[obj]))
    for draw in draws:
        mean = (GMMx[obj,draw], GMMy[obj,draw])
        cov = ((GMMxsig[obj,draw],GMMxycov[obj,draw]),(GMMxycov[obj,draw],GMMysig[obj,draw]))
        hpx, hpy = np.random.multivariate_normal(mean, cov)
        hp_x.append(hpx)
        hp_y.append(hpy)


    
    print(GMMx[obj], GMMy[obj], pi_err[obj])
    

    plt.figure(figsize=(10, 6))
    plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)

    plt.hist2d(hp_x, hp_y, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
    plt.colorbar()
    cmap = cm.get_cmap('viridis')
    rgba = cmap(0)
    ax = plt.axes()
    ax.set_facecolor(rgba)
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.title('Main Sequence', size = 20)
    plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
    plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
#    plt.xlim(xlow, xhigh)
#    plt.ylim(ylow, yhigh)   
    plt.scatter(hp_x, hp_y)
    plt.show()





test = np.load('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/linmix_inputs/linmix_inputs_GMM_2d_4_mTot_delayed/8p5_GMMy.npy')

plt.hist(test)



























