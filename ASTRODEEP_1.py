#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

import csv
import matplotlib.pyplot as plt
import numpy as np

filename_Z  = ['A2744cl_ZPHOT.cat',
                   'A2744PAR_ZPHOT.cat',
                   'M0416cl_ZPHOT.cat',
                   'M0416PAR_ZPHOT.cat',
                   'MACS0717cl_ZPHOT_complete.cat',
                   'MACS0717par_ZPHOT_complete.cat',
                   'MACS1149cl_ZPHOT_complete.cat',
                   'MACS1149par_ZPHOT_complete.cat']

filename_A      = ['A2744cl_A.cat',
                   'A2744PAR_A.cat',
                   'M0416cl_A.cat',
                   'M0416PAR_A.cat',
                   'MACS0717cl_A.cat',
                   'MACS0717par_A.cat',
                   'MACS1149cl_A.cat',
                   'MACS1149par_A.cat']
                  

### GET A COUNT OF Z OBJECTS ###

object_count_z = 0

for i in range(len(filename_Z)):
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_Z[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            object_count_z += 1

object_count_z -= len(filename_Z) # to prevent counting the header in each file



### CREATE ARRAY OF Z DATA FROM CATS FILES ###

data_Z = np.empty([14, object_count_z])
# 0 = filename, 1 = ID, 2 = ZBEST, 3 = object_count (NEVER USE, as HASHED later on)
# 4 onwards are apparent mags filled in later

oc2 = 0 # this object count magically will skip headers..

for i in range(len(filename_Z)):
    
    count = 0
    
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_Z[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            if count == 0:    
                
                # this skips the first line of each file
                count = 1
                
            else:
                
                data_Z[0][oc2] = i
                data_Z[1][oc2] = float(row[0])
                data_Z[2][oc2] = float(row[1])
                data_Z[3][oc2] = oc2
                
                oc2 += 1
                

### REMOVING KNOWN SPURIOUS LINES IN INPUT FILES, redshift of -1 or 99 ###

delete_list = []

for i in range(object_count_z):
    if data_Z[2][i] < 0 or data_Z[2][i] == 99:
        delete_list.append(int(data_Z[3][i]))
        
data_Z = np.delete(data_Z, delete_list, axis=1)
data_Z = np.delete(data_Z, 3, axis=0)


### BACK TO NORMAL

print(len(data_Z[0]))
plt.hist(data_Z[2], bins=10)





### MAGNITUDES ###

### GET A COUNT OF A OBJECTS ###

object_count_A = 0

for i in range(len(filename_Z)):
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_A[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            object_count_A += 1

object_count_A -= len(filename_Z) # to prevent counting the header in each file



### CREATE ARRAY OF MAG DATA FROM CATS FILES ###

data_A = np.empty([16, object_count_A])
# 0 = filename, 1 = ID, 2 onwards are just the colums in the file.


oc3 = 0 # this object count magically will skip headers..

for i in range(len(filename_A)):
    
    count = 0
    
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_A[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            if count == 0:    
                
                # this skips the first line of each file
                count = 1
                
            else:
 
                data_A[0][oc3] = i
                
                for j in range(len(data_A)-1):
                    
                    data_A[j+1][oc3] = float(row[j])

                oc3 += 1
                
print(len(data_A))




### MATCH CAT A AND Z ###

for i in range(len(data_Z[0])):

    for j in range(len(data_A[0])):
        
        if data_Z[0][i] == data_A[0][j] and data_Z[1][i] == data_A[1][j]:
            tester = 4
'''            
            if data_Z[0][i] == 0 or data_Z[0][i] == 1 or data_Z[0][i] == 2 or data_Z[0][i] == 3:
            # these include extra colums X and Y, filters from 5

                for k in range(4, 13):
                    data_Z[k][i] = data_A[k+1][j]

            if data_Z[0][i] == 4 or data_Z[0][i] == 5 or data_Z[0][i] == 6 or data_Z[0][i] == 7:
            # these exclude extra colums X and Y, filters from 3
 
                for k in range(4, 13):
                    data_Z[k][i] = data_A[k-1][j]
'''
print(data_Z)





'''
A2744cl_A.cat
A2744cl_B.cat
A2744cl_C.cat
A2744cl_ZPHOT.cat

A2744PAR_A.cat
A2744PAR_B.cat
A2744PAR_C.cat
A2744PAR_ZPHOT.cat

M0416cl_A.cat
M0416cl_B.cat
M0416cl_C.cat
M0416cl_ZPHOT.cat

M0416PAR_A.cat
M0416PAR_B.cat
M0416PAR_C.cat
M0416PAR_ZPHOT.cat

MACS0717cl_A.cat
MACS0717cl_B.cat
MACS0717cl_ZPHOT_complete.cat

MACS0717par_A.cat
MACS0717par_B.cat
MACS0717par_ZPHOT_complete.cat

MACS1149cl_A.cat
MACS1149cl_B.cat
MACS1149cl_ZPHOT_complete.cat

MACS1149par_A.cat
MACS1149par_B.cat
MACS1149par_ZPHOT_complete.cat
'''