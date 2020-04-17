#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:22:23 2020

@author: lester
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt

param = 'DE'
revision = '105'

samples = 10

# =============================================================================
# OUTPUT - get BEAGLE parameters (<100)
# =============================================================================

fileName = '/Users/lester/Documents/PhD/param_100/fit_{}_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(revision, param)
data_fits = fits.open(fileName)

id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)

data_fits.close()

chi2_fit_arr_allID = [] 

IDs = id_b1
for ID in IDs:
    
    title = '{}-{} ID{} {} samples'.format(param, revision.replace('_', '-'), str(ID), str(samples))
    
    # =============================================================================
    # Get OUTPUT SEDs
    # =============================================================================
    
    data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_{}/{}_BEAGLE.fits.gz'.format(revision, param, ID))
    # print(data_fits.info())
    # print(data_fits['POSTERIOR PDF'].header)
    
    chi2_fit_total = data_fits['POSTERIOR PDF'].data['chi_square']
    
    #needs float64 to provide precision needed for the random.choice weights
    temp_probs = np.float64(data_fits['POSTERIOR PDF'].data['probability'])
    temp_probs = temp_probs/np.sum(temp_probs)
      
    chi2_fit_arr = []
    
    for i in range(samples):
        
        #here's the key line - take weighted samples from the multinest output!
        idx = np.random.choice(len(temp_probs), size=1, p=temp_probs)
        
        # CHI SQUARED
        chi2_fit_arr.append(data_fits['POSTERIOR PDF'].data['chi_square'][idx][0])

    chi2_fit_arr_allID.append(np.average(chi2_fit_arr))

outputDict = {}
outputDict['id']                = id_b1.astype(str)
outputDict['chi2']              = chi2_fit_arr_allID
outputTable = Table(outputDict)
    
#outputTable.write("012_010_chi2.fits", overwrite=True)


plt.hist(outputTable['chi2'])




