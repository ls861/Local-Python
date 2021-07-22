#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 00:06:39 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

directory = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/'
fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
field = fields[0]


for field in fields:
    # BEAGLE OUTPUT CHI SQUARED
    fileName = directory+'chi2/{}_chi2.fits'.format(field)
    data_fits = fits.open(fileName)
    #print(data_fits.info())
    #print(data_fits[1].header)
    id_chi2 = np.asarray(data_fits[1].data['id'], dtype=int)
    chi2 = data_fits[1].data['chi2']
    data_fits.close()
    
    
    print(len(chi2), id_chi2[-1], min(chi2), max(chi2))

