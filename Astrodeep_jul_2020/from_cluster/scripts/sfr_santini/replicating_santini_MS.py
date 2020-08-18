#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:18:30 2020

@author: lester
"""

# Santini approx 1711 sources 1.3 < z < 6.0

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'

AD = np.load(AD_location)
print(AD.dtype.names)
print(len(AD))

# =============================================================================
# H band mag cut of < 27.5
# =============================================================================

AD = AD[AD['H160']<27.5]
print(len(AD))

# =============================================================================
# Photo-z from median of 6 methods
# =============================================================================

zLow = 1.3
zHigh = 2.0

wLim = (zHigh - zLow) / 2.0
zLim = zLow + wLim

AD = AD[abs(AD['ZBEST']-zLim) < wLim]
print(len(AD))

# =============================================================================
# Need sample to be complete above given mass (magnification already included I think)
# =============================================================================

#AD = AD[np.log10(AD['MSTAR']*1e9) > 8.3]
#print(len(AD))

AD = AD[(np.log10(AD['MASTAR_NEB']*1e9)+np.log10(AD['MAGNIF'])) > 8.3]
#AD = AD[(np.log10(AD['MASTAR_NEB']*1e9)) > 8.3]
print(len(AD))

# =============================================================================
# Adding RELFLAG just to be sure...
# =============================================================================

AD = AD[AD['RELFLAG']==1.0]
print('RELFLAG==1', len(AD))

# =============================================================================
# Use high mass cutoff according to Tomczak (between 10.2 up to 10.8 increasing with z)
# =============================================================================

#AD = AD[np.log10(AD['MSTAR']*1e9) < 10.2]
#print(len(AD))

AD = AD[np.log10(AD['MASTAR_NEB']*1e9) < 10.2]
print(len(AD))

# =============================================================================
# get some key columns
# =============================================================================

mass_AD = np.log10(AD['MSTAR']*1e9)
mass_AD_neb = np.log10(AD['MASTAR_NEB']*1e9)
mass_AD_neb_orig = AD['MASTAR_NEB']
sfr_AD = np.log10(AD['SFR'])

field_AD = AD['field']
id_AD = AD['ID']
mag_AD = AD['MAGNIF']


print(field_AD)
print(id_AD)



# =============================================================================
# compare original and BEAGLE IDs for next section
# =============================================================================

#fields = ['A2744_c', 'A2744_p', 'M0416_c', 'M0416_p', 'M0717_c', 'M0717_p', 'M1149_c', 'M1149_p']
#print([fields[1]])
#fields = ['A2744_p', 'M0416_p', 'M0717_p', 'M1149_p']
#fields = fields[0]
#fields = ['A2744_c']
#fields = ['A2744_p']
fields = ['M0416_c']
#fields = ['M0416_p']
#fields = ['M0717_c']
#fields = ['M0717_p']
#fields = ['M1149_c']
#fields = ['M1149_p']


for field in fields:
    fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/data/astrodeep_{}_subset_RF1_001.fits'.format(field)
    data_fits = fits.open(fileName)
    id_input = np.asarray(data_fits[1].data['ID'], dtype=int)
    field_original = np.asarray(data_fits[1].data['field'], dtype=int)
    id_original = np.asarray(data_fits[1].data['ID_original'], dtype=int)
    data_fits.close()
    
    if field == fields[0]:
        print('yes')
        id_input_arr = id_input
        field_original_arr = field_original
        id_original_arr = id_original
        
    else:
        id_input_arr = np.concatenate((id_input_arr, id_input))
        field_original_arr = np.concatenate((field_original_arr, field_original))
        id_original_arr = np.concatenate((id_original_arr, id_original))

    
# =============================================================================
# Santini SFRs (and BEAGLE masses)
# =============================================================================
        
folder_santini = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/scripts/sfr_santini/results/'

mass_BEAGLE_arr = []
sfr_santini_arr = []

for j in range(len(id_AD)):
    
    
    # id_AD and field_AD are AD catalog
    # BEAGLE IDs are needed from id_input_arr, which is same length as id_original and field_original
    # ie match id_AD and field_AD, to id and field ORIGINAL, to then get BEAGLE input id...
    # to do this I need array the size of ORIGINAL / id_input
    

    
    # find BEAGLE input ID, note id_AD and field_AD is new subset
    # field_original_arr only includes selected fields in variable "fields"
    if field_AD[j] in field_original_arr:
        
        arg_field = np.isin(field_original_arr, field_AD[j])
        arg_ID = np.isin(id_original_arr[arg_field], id_AD[j])      

        if sum(arg_ID) != 1:
            print('arg_ID ERROR')
            
#        print(field_original_arr[arg_field][arg_ID], id_original_arr[arg_field][arg_ID], id_input_arr[arg_field][arg_ID])
        
        
        
        '''
    #    arg_field = np.where(field_AD[j]==field_original_arr)[0]
    #    arg_field = np.isin(field_AD[j], field_original_arr)
        arg_field = np.isin(field_original_arr, field_AD[1100])
        print(arg_field)
        print(len(arg_field))  
        print(len(arg_field[arg_field==False]))  
            
        plt.hist(field_original_arr)
        
        
    #    arg_ID = np.where(id_AD[j]==id_original_arr[arg_field])[0]
    #    arg_ID = np.isin(id_original_arr[arg_field], id_AD[j]) 
        
        
        print(len(id_original_arr), len(arg_field))
        print(len(arg_field), sum(arg_field))    
        print(id_AD[j])
        print(id_original_arr[arg_field])    
        print(arg_ID)    
       ''' 
       
   
   
   
        ID_BEAGLE_input = id_input_arr[arg_field][arg_ID]
        
        ID_santini = np.load(folder_santini+'{}_sfr_ID_santini.npy'.format(str(int(field_AD[j]))))
        mass_BEAGLE = np.load(folder_santini+'{}_mass_santini.npy'.format(str(int(field_AD[j]))))
        sfr_santini = np.load(folder_santini+'{}_sfr_santini.npy'.format(str(int(field_AD[j]))))
        sfr_err_santini = np.load(folder_santini+'{}_sfr_err_santini.npy'.format(str(int(field_AD[j]))))
    

        idx = np.isin(ID_santini, ID_BEAGLE_input)
        
        if sum(idx) == 0:
            mass_BEAGLE_arr.append(-11.0)
            sfr_santini_arr.append(-11.0)
            print('1 ERROR')
        elif sum(idx) == 1:
            mass_BEAGLE_arr.append(mass_BEAGLE[idx][0] - np.log10(mag_AD[j]))
            sfr_santini_arr.append(sfr_santini[idx][0] - np.log10(mag_AD[j]))
#            sfr_santini_arr.append(sfr_santini[idx][0] - 0)
        else:
            print('2 ERROR')
            
    else:
        mass_BEAGLE_arr.append(-10.0)
        sfr_santini_arr.append(-10.0)
        print('0 ERROR')       

        
##    idx = np.where(ID_BEAGLE_input==ID_santini)[0]
#    idx = np.isin(ID_santini, ID_BEAGLE_input)
#    if sum(idx) == 0:
#        mass_BEAGLE_arr.append(-10.0)
#        sfr_santini_arr.append(-10.0)
#        print('0 ERROR')
#    elif sum(idx) == 1:
#        mass_BEAGLE_arr.append(mass_BEAGLE[idx][0] - np.log10(mag_AD[j]))
#        sfr_santini_arr.append(sfr_santini[idx][0] - np.log10(mag_AD[j]))
#    else:
#        print('ERROR')




print(len(id_AD), len(mass_BEAGLE_arr), len(ID_santini))


'''
Summary so far...
AD: mass_AD, mass_AD_neb, sfr_AD, SFR_NEB
BEAGLE: mass_BEAGLE_arr
Santini: sfr_santini_arr
'''


field_AD = np.array(field_AD)
id_AD = np.array(id_AD)
mag_AD = np.array(mag_AD)

mass_AD = np.array(mass_AD)
mass_AD_neb = np.array(mass_AD_neb)
mass_AD_neb_orig = np.array(mass_AD_neb_orig)
mass_BEAGLE_arr = np.array(mass_BEAGLE_arr)
sfr_AD = np.array(sfr_AD)
sfr_santini_arr = np.array(sfr_santini_arr)

idx = (mass_BEAGLE_arr>0.0) # these are the ones which santini ID didn't match AD ID, because of initial filters 

print(len(mass_AD), len(mass_AD_neb), len(mass_BEAGLE_arr), len(sfr_AD), len(sfr_santini_arr))

#plt.hist(sfr_santini_arr)
#plt.show()
#plt.hist(mass_AD_neb)
#plt.show()

field_AD = field_AD[idx]
id_AD = id_AD[idx]
mag_AD = mag_AD[idx]

mass_AD = mass_AD[idx]
mass_AD_neb = mass_AD_neb[idx]
mass_AD_neb_orig = mass_AD_neb_orig[idx]
mass_BEAGLE_arr = mass_BEAGLE_arr[idx]
sfr_AD = sfr_AD[idx]
sfr_santini_arr = sfr_santini_arr[idx]

#plt.hist(sfr_santini_arr)
#plt.show()
#plt.hist(mass_AD_neb)
#plt.show()

print(min(sfr_santini_arr))

# =============================================================================
# some plots
# =============================================================================

#plt.scatter(mass_AD, sfr_AD)
#plt.show()
#
#plt.scatter(mass_BEAGLE_arr, sfr_AD)
#plt.show()
#
#plt.scatter(mass_AD, sfr_santini_arr)
#plt.show()
#
#plt.scatter(mass_BEAGLE_arr, sfr_santini_arr)
#plt.show()
#
#plt.scatter(mass_AD, mass_BEAGLE_arr)
#plt.show()
#
#plt.scatter(sfr_AD, sfr_santini_arr)
#plt.show()

# =============================================================================
# 2 sigma clipping to remove outliers
# =============================================================================

#mass = mass_BEAGLE_arr
#sfr = sfr_santini_arr

#mass = mass_AD
#sfr = sfr_AD

mass = mass_AD_neb
sfr = sfr_santini_arr

#mass = mass_AD_neb
#sfr = sfr_AD

outliers = 1

while outliers > 0:
    
    fit = np.polyfit(mass, sfr, 1)
    print(fit)

    sfr_residuals= sfr - (fit[0]*mass + fit[1])
    
    sigma = np.std(sfr_residuals)
    
    idx = (abs(sfr_residuals)<2*sigma)
    
    outliers = len(mass) - sum(idx)

    mass = mass[idx]
    sfr = sfr[idx]


# =============================================================================
# Santini
# =============================================================================

# logSFR = alpha log(M / M_9.7) + beta

z_med_san = np.array((1.65, 2.5, 3.5, 4.5, 5.5))
A_san = np.array((1.04, 1.16, 1.02, 0.94, 0.92))
A_err_san = np.array((0.03, 0.03, 0.04, 0.06, 0.15))
B_san = np.array((1.01, 1.22, 1.37, 1.37, 1.99))
B_err_san = np.array((0.04, 0.03, 0.03, 0.05, 0.13))

alpha_san = B_san - 9.7*A_san
beta_san = A_san

alpha_err_san = (B_err_san**2 + (9.7*A_err_san)**2) ** 0.5
beta_err_san = A_err_san

x = np.linspace(min(mass), max(mass))

plt.scatter(mass, sfr)
plt.plot(x, fit[0]*x + fit[1], color='k', label='lester')
plt.plot(x, beta_san[0]*x + alpha_san[0], color='r', label='santini')
plt.title(field.replace('_',''))
plt.legend()
plt.show()

print(fit[0], fit[1])



#print(id_AD[66])
#print(mass_AD_neb[66])
#print(mass_AD_neb_orig[66])
#
#print(mass_AD[66])
#print((10**mass_AD[66])/1e9)
#print(field_AD[66])
#
#print(np.log10(4.075*1e9))
#
#
#print(mag_AD[66])








