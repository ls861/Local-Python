#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:22:23 2020

@author: lester
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
#import matplotlib.pyplot as plt

params = ['0A2744C', '1A2744P', '2M0416C', '3M0416P']
params = ['4M0717C', '5M0717P', '6M1149C', '7M1149P']
revision = '001'

for param in params:
    
    # =============================================================================
    # BEAGLE OUTPUT
    # =============================================================================
    
    fileName = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(param, revision)
    
    data_fits = fits.open(fileName)
    
    id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    
    data_fits.close()
    
    chi2_fit_arr_allID = [] 
    
    IDs = id_b1
    
    for ID in IDs:
#        print(str(ID) + '/' + str(len(IDs)))
        
#        title = '{}-{} ID{} {} samples'.format(param, revision.replace('_', '-'), str(ID), str(samples))
        
        # =============================================================================
        # Get min ch2
        # =============================================================================
        
        data_fits = fits.open('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_{}/{}_BEAGLE.fits.gz'.format(param, revision, ID))
        # print(data_fits.info())
        # print(data_fits['POSTERIOR PDF'].header)
        
        chi2_fit_total = data_fits['POSTERIOR PDF'].data['chi_square']
        data_fits.close()
        
        chi2_fit_arr_allID.append(np.amin(chi2_fit_total))
    
    outputDict = {}
    outputDict['id']                = id_b1.astype(str)
    outputDict['chi2']              = chi2_fit_arr_allID
    outputTable = Table(outputDict)
        
    outputTable.write("/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_run/{}_chi2.fits".format(param), overwrite=True)
    
    #plt.hist(outputTable['chi2'], bins=100)
    #plt.ylim(0, 10)
    #plt.show()
    #
    #plt.hist(outputTable['chi2'], bins=100, range=(0, 200))
    #plt.show()
    #
    #print(len(outputTable['chi2']))
    #print(len(outputTable['chi2'][outputTable['chi2']<=25]))
    #print(len(outputTable['chi2'][outputTable['chi2']>25]))
    
    




