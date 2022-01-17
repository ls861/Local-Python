#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 11:55:04 2021

@author: lester
"""

import numpy as np
import pickle
from astropy.io import fits
from redshift_vs_redshift import redshift_vs_redshift


figuresize = 7
fontsize_legend = 12
fontsize_axes = 20
save = False

# AD catalogue
AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/new_beagle/data/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)

with open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/data.p', 'rb') as f:
    ADx = pickle.load(f, encoding='latin1')

folder = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/BEAGLE_dust_install/analysis/kelly_input/'

fileNames = ['scenario_34_clusters_z1p25-6p0.fits']
s34 = fits.open(folder+fileNames[0])[1].data


redshift_vs_redshift(AD, ADx, s34, figuresize, fontsize_legend, fontsize_axes, save)






