#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 01:22:36 2020

@author: lester
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def mag(flux):
    '''
    Takes flux in uJy and returns ABmag
    '''
    ABmag = 23.9 - 2.5*np.log10(flux)
    return ABmag



def flux(mag):
    '''
    Takes abs(mag) and returns flux in uJy 
    '''
    flux =10 ** ((23.9 - np.abs(mag)) / 2.5)
    return flux


# dtype=[('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]

D = np.load('astrodeep_rawfile.npy')

D0 = D[D['field']==1]

#print(len(D))
#print(len(D0))

# TEST OF BASIC ONE

print(D0['b_H160'][0])
print(D0['H160'][0])
test_mag = mag(D0['b_H160'][0])
print(test_mag)
test_flux = flux(D0['H160'][0])
print(test_flux)

# NEGATIVE FLUXES

print('FLUX and ERROR')
D0mf = D0[D0['b_H160']<0] # field 0, minus flux
print(D0mf['b_H160'])
print(D0mf['b_errH160'])


upper_D0mf = D0mf['b_H160'] + D0mf['b_errH160']

print('UPPER FLUX')
print(upper_D0mf)

print('FLUX and ERROR - CALCULATED FROM MAGS')
print(flux(D0mf['H160']))
print(flux(D0mf['errH160']))

print('OTHER')
print(D0mf['ID'])
print(D0mf['RELFLAG'])
print(D0mf['a_ID'])
print(D0mf['b_ID'])
print(D0mf['H160'])
print(D0mf['errH160'])


# PLOT

plt.scatter(abs(D0mf['b_H160']), flux(D0mf['H160']), marker='x')
plt.scatter(abs(upper_D0mf), flux(D0mf['H160']))
plt.plot((0, 0.1), (0, 0.1))
plt.xlim(0, 0.02)
plt.ylim(0, 0.02)

plt.show()













# CALCULATED MAGS

#print(' ')
#print('CALCULATED MAGS')
#D0mm = mag(D0mf['b_H160'])
#upper_D0mm = mag(upper_D0mf)
#
#print(D0mm)
#print(upper_D0mm)


# MAGS FROM CAT

#print('   ')
#D0mm = D0[D0['H160']<0] # field 0, minus mag
#print(D0mm['b_H160'])
#print(D0mm['b_errH160'])
#print(D0mm['H160'])
#print(D0mm['errH160'])





# =============================================================================
# theory test
# =============================================================================


test = np.load('astrodeep_rawfile.npy')

test = test[test['field'] < 4]
test = test[test['ID'] < 20000]
test = test[test['b_H160'] < 0]

print(test['field'])
print(test['b_H160'])
print(test['H160'])

print(len(test))

test = np.load('astrodeep_rawfile.npy')

test = test[test['field'] < 4]
test = test[test['H160'] > 90]

print(len(test))












