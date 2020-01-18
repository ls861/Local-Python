#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 18:45:55 2019

@author: lester
"""


import numpy as np

flux = 343543543 # uJ
apparent_mag = -2.5 * np.log10(  (flux * (10**-6))  /  3631) # calculated apparent mags from flux

print(flux)
print(apparent_mag)


magnification = 7


int_flux = flux / magnification

#1 (amend flux)
int_apparent_mag_1 = -2.5 * np.log10(  (int_flux * (10**-6))  /  3631) # calculated apparent mags from flux

#2 (amend mag)
int_apparent_mag_2 = apparent_mag + (2.5 * np.log10(magnification))



print(int_flux)
print(int_apparent_mag_1)
print(int_apparent_mag_2)