#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:36:19 2020

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits



beagleData = fits.open('/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_006/astrodeep_002/pyp-beagle/data/BEAGLE_summary_catalogue.fits')



#needs float64 to provide precision needed for the random.choice weights
temp_probs = np.float64(beagleData['POSTERIOR PDF'].data['probability'])
 temp_probs = temp_probs/np.sum(temp_probs)
 mass = []
 sfr = []
 for j in range(1000):
    #here’s the key line - take weighted samples from the multinest output!
    idx = np.random.choice(len(temp_probs), p=temp_probs)
    mass.append(np.log10(beagleData['GALAXY PROPERTIES'].data['M_tot'][idx]))
    sfr.append(np.log10(beagleData['STAR FORMATION'].data['sfr'][idx]))














