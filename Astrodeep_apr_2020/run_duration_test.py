#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 18:09:05 2020

@author: lester
"""

from astropy.io import fits

param = 'DE'
revision = '106'
ID = 1

data_fits = fits.open('/Users/lester/Documents/PhD/param_100/fit_{}_{}/{}_BEAGLE.fits.gz'.format(revision, param, ID))
#print(data_fits.info())
print(data_fits[0].header)
data_fits.close()

    







































