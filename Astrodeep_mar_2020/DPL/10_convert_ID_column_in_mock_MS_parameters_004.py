#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 10:28:32 2020

@author: lester
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np


# =============================================================================
# get INPUT params (len 100)
# =============================================================================

outputDict = {}

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_mar_2020/DPL/mock_MS_parameters_004.fits'
data_fits = fits.open(fileName)
print(data_fits[1].header)

ID                              = data_fits[1].data['id']
outputDict['id']                = ID.astype(str)
outputDict['tau']               = data_fits[1].data['tau']
outputDict['dpl_alpha']         = data_fits[1].data['dpl_alpha']
outputDict['dpl_beta']          = data_fits[1].data['dpl_beta']
outputDict['metallicity']       = data_fits[1].data['metallicity']
outputDict['mass']              = data_fits[1].data['mass']

outputDict['nebular_logU']      = data_fits[1].data['nebular_logU']
outputDict['nebular_xi']        = data_fits[1].data['nebular_xi']

outputDict['tauV_eff']          = data_fits[1].data['tauV_eff']

outputDict['A']                 = data_fits[1].data['A']
outputDict['sfr']               = data_fits[1].data['sfr']
outputDict['closest_ind']       = data_fits[1].data['closest_ind']

outputDict['ssfr']              = np.log10(data_fits[1].data['sfr']) - data_fits[1].data['mass']

data_fits.close()

outputTable = Table(outputDict)
    
outputTable.write("mock_MS_parameters_005.fits", overwrite=True)









