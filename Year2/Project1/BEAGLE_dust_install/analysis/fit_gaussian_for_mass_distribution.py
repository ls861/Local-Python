#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 17:53:02 2021

@author: lester
"""

import numpy as np
from scipy.stats import norm
import pickle
import matplotlib.pyplot as plt


scenarioA = '27'
fields = ['clusters']
z_bins = ['z1p25-6p0','z1p25-2p0','z2p0-3p0','z3p0-4p0','z4p0-5p0']
z_bins = ['z5p0-6p0']
#z_bins = ['z1p25-6p0']

for field in fields:
    for z_bin in z_bins:
        
        data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/scenario_{}_{}_{}.p'.format(scenarioA, field, z_bin),'r'))
        
        data = data['mass_BEAGLE_stellar']
        mean,std=norm.fit(data)
                
                
                
        plt.hist(data, bins=30, normed=True)
        plt.title('{} {} {}, mean:{:.2f} std:{:.2f}'.format(scenarioA, field, z_bin, mean, std))
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        y = norm.pdf(x, mean, std)
        plt.plot(x, y)
        plt.show()
        
        print(mean, std)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        