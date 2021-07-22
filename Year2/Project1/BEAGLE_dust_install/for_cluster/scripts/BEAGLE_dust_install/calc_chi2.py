#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 16:22:23 2020
@author: lester
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
import os
#import matplotlib.pyplot as plt

subfolder = 'BEAGLE_dust_install'
params = ['0', '1', '2', '3', '4', '5', '6', '7']
revision = '001'

for param in params:

    # =============================================================================
    # BEAGLE OUTPUT
    # =============================================================================

    fileName = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/{}/fit_{}/pyp-beagle/data/BEAGLE_summary_catalogue.fits'.format(subfolder, param, revision)
    data_fits = fits.open(fileName)
    id_b1 = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    data_fits.close()
    min_chi2_arr = []
    IDs = id_b1.astype(str)

    for ID in IDs:
#        print(str(ID) + '/' + str(len(IDs)))
#        title = '{}-{} ID{} {} samples'.format(param, revision.replace('_', '-'), str(ID), str(samples))

        # =============================================================================
        # Get min ch2
        # =============================================================================
        data_fits = fits.open('/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/{}/fit_{}/{}_BEAGLE.fits.gz'.format(subfolder, param, revision, ID))
        # print(data_fits.info())
        # print(data_fits['POSTERIOR PDF'].header)
        chi2 = data_fits['POSTERIOR PDF'].data['chi_square']
        data_fits.close()
        min_chi2_arr.append(np.amin(chi2))

    outputDict = {}
    outputDict['id']                = IDs
    outputDict['chi2']              = min_chi2_arr
    outputTable = Table(outputDict)

    os.system("mkdir ./chi2")
    outputTable.write("./chi2/{}_chi2.fits".format(param), overwrite=True)

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
