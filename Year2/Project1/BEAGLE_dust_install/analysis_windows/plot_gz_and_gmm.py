#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 16:51:46 2021

@author: lester
"""


import numpy as np
import pickle
import matplotlib.pyplot as plt



scenarioA = '29'
field = 'clusters'
z_bin = 'z1p25-2p0'
z_lower = 1.25
z_upper = 2.0

#WINDOWS
with open('/Users/LSand/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis_windows/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin), 'rb') as f:
    s = pickle.load(f, encoding='latin1')

field_AD = 0
id_AD = 2359

z_med_hp = (z_lower+z_upper)/2.0
z_med_hp_gap = (z_lower+z_upper)/2.0 - z_lower

x_hp = np.array([])
y_hp = np.array([])
z_hp = np.array([])

n_hp = 3000 # number of samples to take from GMM in total

for i in range(len(s['id_AD'])):
    
    
#    print(i, s['field_AD'][i], field_AD, s['id_AD'][i], id_AD)
    
#    if s['field_AD'][i] == field_AD and s['id_AD'][i] == id_AD:
#    if i in [212, 259, 171, 298, 382,   7,  87, 389, 384,  45]: #mass
#    if i in [175,  58, 249, 254, 309, 171, 150,  60, 100, 382]: #sfr
#    if i in [100, 116, 298,  51,  64,  99, 350, 112, 171, 347]: #redshift
    if i in [350]: #redshift        
        x_temp = np.array([])
        y_temp = np.array([])
        z_temp = np.array([])
    
        for G in range(3):
            
            mean = np.array([s['x_GMM_3d'][i,G],s['y_GMM_3d'][i,G],s['z_GMM_3d'][i,G]])
            cov = np.array([[np.power(s['xsig_GMM_3d'][i,G],2), s['xycov_GMM_3d'][i,G], s['xzcov_GMM_3d'][i,G]],[s['xycov_GMM_3d'][i,G], np.power(s['ysig_GMM_3d'][i,G],2), s['yzcov_GMM_3d'][i,G]],[s['xzcov_GMM_3d'][i,G], s['yzcov_GMM_3d'][i,G], np.power(s['zsig_GMM_3d'][i,G],2)]])
    
            xyz = np.random.multivariate_normal(mean, cov, size=int(n_hp*s['amp_GMM_3d'][i,G]))
    
            x_hp = np.concatenate((x_hp,xyz[:,0]))
            y_hp = np.concatenate((y_hp,xyz[:,1]))
            z_hp = np.concatenate((z_hp,xyz[:,2]))
    
            x_temp = np.concatenate((x_temp,xyz[:,0]))
            y_temp = np.concatenate((y_temp,xyz[:,1]))
            z_temp = np.concatenate((z_temp,xyz[:,2]))
            
        plt.scatter(x_temp,y_temp,c=z_temp)
        plt.title('{} {} {} {}'.format(i, int(s['field_AD'][i]), int(s['id_AD'][i]), int(s['id_BEAGLE'][i])))
        plt.colorbar()
        plt.show()

        plt.scatter(x_temp,z_temp,c=y_temp)
        plt.title('{} {} {}'.format(i, int(s['field_AD'][i]), int(s['id_AD'][i])))
        plt.colorbar()
        plt.show()
        
        plt.scatter(y_temp,z_temp,c=x_temp)
        plt.title('{} {} {}'.format(i, int(s['field_AD'][i]), int(s['id_AD'][i])))
        plt.colorbar()
        plt.show()



















#%%
i = 314
piBeagle = s['amp_GMM_3d']
maxind = np.argmax(piBeagle[i])
rdmind = np.random.choice(3, p=piBeagle[i])
print(piBeagle[i], maxind, rdmind)



