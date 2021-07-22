#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 18:33:49 2020

@author: lester
"""

import os
import numpy as np
import matplotlib.pyplot as plt

folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'


fileList = os.listdir(folder)
for f, file in enumerate(fileList):
    if '.fits.gz' in file:
        
        ID = file.replace('_BEAGLE.fits.gz','')
    
        sfr = np.load(folder+'{}_log_instant_sfr.npy'.format(ID))
    
    
        print(sfr)
    
    
        plt.hist(sfr, bins=30)
        plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        




