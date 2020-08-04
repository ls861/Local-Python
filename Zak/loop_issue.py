#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 08:38:11 2020

@author: lester
"""



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

IRAS_03254_Spec = pd.read_csv('i03254_YSO.csv')
L1014_BG_Spec = pd.read_csv('L1014_BG_Star.csv')
L1489_Spec = pd.read_csv('L1489_YSO.csv')
SVS4_5_Spec = pd.read_csv('SVS4-5_YSO.csv')

#IRAS_03254_Phot = pd.read_csv('i03254_Photo.csv')
#L1014_BG_Phot = pd.read_csv('L1014_Photo.csv')
#L1489_Phot = pd.read_csv('L1489IRS_Photo.csv')
#SVS4_5_Phot = pd.read_csv('SVS4-5_Photo.csv')

IRAS_03254_355 = 0 
L1014_BG_355 = 0
L1489_355 = 0
SVS4_5_355 = 0

objectNames = ['IRAS_03254','L1014_BG','L1489','SVS4_5']
spec_list = [IRAS_03254_Spec,L1014_BG_Spec,L1489_Spec,SVS4_5_Spec]
spec_df = pd.concat(spec_list,keys=['IRAS_03254','L1014_BG','L1489','SVS4_5'])

#%%

#spec_df.loc['IRAS_03254']

spec_df = spec_df[["Wavelength(um)",'Flux(Jy)']]
spec_df['Flux(mJy)'] = spec_df['Flux(Jy)']*1000
spec_df = spec_df[["Wavelength(um)",'Flux(mJy)']]
print(spec_df)



#%%

spec_df.loc['SVS4_5']


test = spec_df.loc[objectNames[3]].drop_duplicates(subset="Wavelength(um)")
print(test)


test = [0, 0, 0, 0]
 
  
  
test[0] = spec_df.loc[objectNames[3]].drop_duplicates(subset="Wavelength(um)")
print(test[0])
  
  
for i in range(len(objectNames)):
  test[i] = spec_df.loc[objectNames[3]].drop_duplicates(subset="Wavelength(um)")



print(test[0])
  

spec_df.loc[objectNames[3]]


i = 0

lester                      = spec_df.loc[objectNames[i]].drop_duplicates(subset="Wavelength(um)")
spec_df.loc[objectNames[i]] = spec_df.loc[objectNames[i]].drop_duplicates(subset="Wavelength(um)")

print(lester)
print(spec_df.loc[objectNames[i]])



spec_df.loc['SVS4_5']
spec_df = spec_df.drop_duplicates(subset="Wavelength(um)")
spec_df.loc['SVS4_5']

spec_df['SVS4_5'] = 

print(type(spec_df.loc['SVS4_5']))




for index in range(len(objectNames)):  
    #  This for loop has organised the spectral data in a format which is usable in the ETC WITHOUT SCALING TO 44.6mJy.
    spec_df.loc[objectNames[index]] = spec_df.loc[objectNames[index]].drop_duplicates(subset="Wavelength(um)")
    
    #spec_df.loc[ObjectName] = spec_df.loc[ObjectName].reset_index(drop=True)
    #spec_list[index].to_csv(f"{ObjectName}.dat", header=False, index=False)


spec_df.loc[objectNames[2]]


pri



























