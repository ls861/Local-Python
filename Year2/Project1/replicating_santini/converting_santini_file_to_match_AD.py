#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:35:13 2020

@author: lester
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

real_SAN_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/replicating_santini/final_cat_Santini17.dat'
columns_to_keep = ['#ID', 'z', 'lgM', 'lgSFR', 'magnif']
df = pd.read_csv(real_SAN_location, sep="\s+", usecols=columns_to_keep)

#plt.hist(df['#ID'].values, bins=100)
#plt.show()

AD_location = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/astrodeep/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)

id_SANTINI = []
mass_SANTINI = []
sfr_SANTINI = []
redshift_SANTINI = []
mag_SANTINI = []


for i in range(len(AD)):
#for i in range(1):
    
    if AD['field'][i]==0:
        id_adjustment = 0.0
    elif AD['field'][i]==2:
        id_adjustment = 3e4
    elif AD['field'][i]==4:
        id_adjustment = 4e4
    elif AD['field'][i]==6:
        id_adjustment = 7e4
    else:
        id_adjustment = -3e5
        
    temp_id_santini = AD['ID'][i] + id_adjustment    
    temp_idx = np.isin(df['#ID'].values, temp_id_santini)
                         
    if sum(temp_idx) == 0:
        
        id_SANTINI.append(-103.0)
        mass_SANTINI.append(-103.0)
        sfr_SANTINI.append(-103.0)
        redshift_SANTINI.append(-103.0)
        mag_SANTINI.append(-103.0)
        
    elif sum(temp_idx) == 1:

        id_SANTINI.append(df['#ID'].values[temp_idx])
        mass_SANTINI.append(df['lgM'].values[temp_idx])
        sfr_SANTINI.append(df['lgSFR'].values[temp_idx])
        redshift_SANTINI.append(df['z'].values[temp_idx])
        mag_SANTINI.append(df['magnif'].values[temp_idx])    

    else:
        print('LEARN TO CODE')
                          
       
id_SANTINI = np.array(id_SANTINI)
mass_SANTINI = np.array(mass_SANTINI)
sfr_SANTINI = np.array(sfr_SANTINI)
redshift_SANTINI = np.array(redshift_SANTINI)
mag_SANTINI = np.array(mag_SANTINI)

#np.save('id_SANTINI',id_SANTINI)
#np.save('mass_SANTINI',mass_SANTINI)
#np.save('sfr_SANTINI',sfr_SANTINI)
#np.save('redshift_SANTINI',redshift_SANTINI)
#np.save('mag_SANTINI',mag_SANTINI)

















                   
                    

