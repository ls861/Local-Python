#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 12:56:05 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# [('field', '<f8'), ('ID', '<f8'), ('ZBEST', '<f8'), ('ZBEST_SIQR', '<f8'), ('MAGNIF', '<f8'), ('ZSPECFLAG', '<f8'), ('chi2', '<f8'), ('MSTAR', '<f8'), ('MASTAR_MIN', '<f8'), ('MSTAR_MAX', '<f8'), ('SFR', '<f8'), ('SFR_MIN', '<f8'), ('SFR_MAX', '<f8'), ('chi2_NEB', '<f8'), ('MASTAR_NEB', '<f8'), ('MASS_MIN_NEB', '<f8'), ('MASS_MAX_NEB', '<f8'), ('SFR_NEB', '<f8'), ('SFR_MIN_NEB', '<f8'), ('SFR_MAX_NEB', '<f8'), ('RELFLAG', '<f8'), ('a_ID', '<f8'), ('RA', '<f8'), ('DEC', '<f8'), ('B435', '<f8'), ('V606', '<f8'), ('I814', '<f8'), ('Y105', '<f8'), ('J125', '<f8'), ('JH140', '<f8'), ('H160', '<f8'), ('Ks', '<f8'), ('CH1', '<f8'), ('CH2', '<f8'), ('errB435', '<f8'), ('errV606', '<f8'), ('errI814', '<f8'), ('errY105', '<f8'), ('errJ125', '<f8'), ('errJH140', '<f8'), ('errH160', '<f8'), ('errKs', '<f8'), ('errCH1', '<f8'), ('errCH2', '<f8'), ('b_ID', '<f8'), ('b_B435', '<f8'), ('b_V606', '<f8'), ('b_I814', '<f8'), ('b_Y105', '<f8'), ('b_J125', '<f8'), ('b_JH140', '<f8'), ('b_H160', '<f8'), ('b_Ks', '<f8'), ('b_CH1', '<f8'), ('b_CH2', '<f8'), ('b_errB435', '<f8'), ('b_errV606', '<f8'), ('b_errI814', '<f8'), ('b_errY105', '<f8'), ('b_errJ125', '<f8'), ('b_errJH140', '<f8'), ('b_errH160', '<f8'), ('b_errKs', '<f8'), ('b_errCH1', '<f8'), ('b_errCH2', '<f8')]
# =============================================================================

'''
catalogs_b        = [A2744_c_b,
                     A2744_p_b,
                     M0416_c_b,
                     M0416_p_b,
                     M0717_c_b,
                     M0717_p_b,
                     M1149_c_b,
                     M1149_p_b]
'''
# take random sample from astrodeep, uniformly sampled in redshift and Hband mag, maybe 200? then can plot the times - POTENTIALLY TAKE 1 PARALLEL


D = np.load("astrodeep_rawfile.npy")
D1 = np.sort(D[D['field'] == 1], order='ZBEST')

ID_new = []

for i in range(10):
    D1_z_subset = np.sort(D1[(D1['ZBEST'] > i) & (D1['ZBEST'] <= (i+1))], order='b_H160')
    j = 20 # number per redshift bin to sample
    if len(D1_z_subset) < j:
        j = len(D1_z_subset)
    ID_new = ID_new + list((D1_z_subset[::int(len(D1_z_subset)/j)]['ID']))
    
idx = []

for i in range(len(D1)):
    if D1['ID'][i] in ID_new:
        idx.append(i)

D1_new = D1[idx]

plt.hist(D1_new['ZBEST'])
plt.show()

plt.hist(D1_new['b_H160'], bins=100, range=(0, 25))
plt.show()

print(len(D1_new))


plt.scatter(D1_new['ZBEST'], D1_new['b_H160'])
plt.xlim(0, 10)
plt.ylim(0, 0.5)
plt.show()


# =============================================================================
# filters = ['HST_ACS_WFC_F435W', 'HST_ACS_WFC_F606W', 'HST_ACS_WFC_F814W', 'HST_WFC3_IR_F105W', 'HST_WFC3_IR_F125W', 'HST_WFC3_IR_F140W', 'HST_WFC3_IR_F160W', 'Paranal_HAWKI_Ks', 'Spitzer_IRAC_I1', 'Spitzer_IRAC_I2']
# =============================================================================

header = ['b_B435', 'b_errB435', 'b_V606', 'b_errV606', 'b_I814', 'b_errI814', 'b_Y105', 'b_errY105', 'b_J125', 'b_errJ125', 'b_JH140', 'b_errJH140', 'b_H160', 'b_errH160', 'b_Ks', 'b_errKs', 'b_CH1', 'b_errCH1', 'b_CH2', 'b_errCH2', 'ZBEST', 'field', 'ID']

header_string = '#ID b_B435 b_errB435 b_V606 b_errV606 b_I814 b_errI814 b_Y105 b_errY105 b_J125 b_errJ125 b_JH140 b_errJH140 b_H160 b_errH160 b_Ks b_errKs b_CH1 b_errCH1 b_CH2 b_errCH2 ZBEST field ID_original\n'


#f= open("astrodeep_A2744_p_subset_001.txt","w+")
#f.write(header_string)
#
#for i in range(len(D1_new)):
#    row = str(i+1) # ID from 1 to len(D)
#    
#    for j in range(len(header)):
#        row = row + ' ' + str(D1_new[header[j]][i])
#        
#    row = row + '\n'    
#    f.write(row)
#
#f.close()

#np.save("D1_new", D1_new)











