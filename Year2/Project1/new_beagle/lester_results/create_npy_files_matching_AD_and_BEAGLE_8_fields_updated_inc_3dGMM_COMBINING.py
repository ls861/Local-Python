#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:31:27 2020

@author: lester
"""

import numpy as np
import pickle

#subfolder = 'delayed_uniform_logtau'
#subfolder = 'danger_delayed_uniform_logtau'
#subfolder = 'danger_constant'
subfolder = 'safety'

names = ['id_BEAGLE', 'mass_BEAGLE_tot', 'mass_BEAGLE_stellar', 'sfr_BEAGLE_instant', 'redshift_BEAGLE', 'tau_BEAGLE', 'tauv_BEAGLE', 'msa_BEAGLE', 'metallicity_BEAGLE', 'min_chi2_BEAGLE', 'id_GMM_2d', 'x_GMM_2d', 'y_GMM_2d', 'xsig_GMM_2d', 'ysig_GMM_2d', 'xycov_GMM_2d', 'amp_GMM_2d', 'id_GMM_3d', 'x_GMM_3d', 'y_GMM_3d', 'z_GMM_3d', 'xsig_GMM_3d', 'ysig_GMM_3d', 'zsig_GMM_3d', 'xycov_GMM_3d', 'xzcov_GMM_3d', 'yzcov_GMM_3d', 'amp_GMM_3d']

dic = {}

for name in names:

    dic[name] = np.load('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files/{}_{}_0.npy'.format(subfolder, name))

    for i in range(8-1):
        if len(dic[name].shape) == 1:
            dic[name] = np.hstack((dic[name], np.load('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files/{}_{}_{}.npy'.format(subfolder, name, i+1))))
        else:
            dic[name] = np.vstack((dic[name], np.load('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/lester_results/npy_files/{}_{}_{}.npy'.format(subfolder, name, i+1))))

pickle.dump(dic, open('./{}/{}_pickle.p'.format(subfolder, subfolder),'w')) 




