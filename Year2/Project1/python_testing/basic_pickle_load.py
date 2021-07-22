#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 17:59:40 2020

@author: lester
"""

import pickle 






data = pickle.load(open('/Users/lester/Documents/BEAGLE_heatplot_samples/0_23_BEAGLE_samples.p','r'))
            
            
#data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/report_BEAGLE_instant_SFRs/{}_instant_sfr_medians.p'.format(int(field)),'r'))       
            
            
print(data)
         
import pickle 
for i in range(8):
    data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/GMM/kelly_2d_GMM_inputs/{}_mStar_delayed_GMM_2d_fits.p'.format(i),'r'))

    print(len(data['id']))


print(data.keys())
print(data['id'])
print(data['nGauss'])
