#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import math

np.set_printoptions(threshold=sys.maxsize)



filename_Z      = ['A2744cl_ZPHOT.cat',
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
            if row[0] == '#':
                pass
            elif int(row[0]) > 100000:
                pass
            else:
                object_count_z += 1 # 29373


### GET A COUNT OF A OBJECTS ###

object_count_A = 0

for i in range(len(filename_Z)):
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_A[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            if row[0] == '#' or row[0] == '#ID':
                pass
            elif int(row[0]) > 100000:
                pass
            else:
                object_count_A += 1 # 29373

### BASIC CHECK TO CONFIRM NO MAJOR ERRORS ###

if object_count_A == object_count_z:
    print('The number of objects in A and Z are consistent: ' + str(object_count_A))
else:
    print("You've ballsed up somewhere Lester")



### CREATE ARRAY OF Z DATA FROM CATS FILES ###

data_Z = np.empty([15, object_count_z]) 

# note initialises with random entries...
'''
0   filename
1   ID
2   ZBEST
3   object_count (NEVER USE, as HASHED later on)
4   ID from mag to check consistency
5   MAG_B435
6   MAG_V606 
7   MAG_I814 
8   MAG_Y105 
9   MAG_J125 
10  MAG_JH140 
11  MAG_H160 
12  MAG_Ks 
13  MAG_IRAC1 
14  MAG_IRAC2

'''

oc2 = 0 # this object count magically will skip headers..

for i in range(len(filename_Z)):
    
    count = 0
    
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_Z[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            if count == 0:    
                
                # this skips the first line of each file
                count = 1
                
            elif int(row[0]) > 100000:
                pass
                
            else:
                
                data_Z[0][oc2] = i
                data_Z[1][oc2] = float(row[0])
                data_Z[2][oc2] = float(row[1])
                data_Z[3][oc2] = oc2
                
                oc2 += 1
                

### INCORPORATE MAG DATA FROM CATS FILES ###

# added to columns 4 to 13 (ID + 10 filters)

oc3 = 0 # this object count magically will skip headers..

for i in range(len(filename_A)):
    
    count = 0
    
    with open('/Users/lester/BEAGLE/ASTRODEEP/' + filename_A[i] ) as csvfile:
        fits = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in fits:
            
            if count == 0:    
                
                # this skips the first line of each file
                count = 1
                
            elif int(row[0]) > 100000:
                pass 
            
            else:
                
                data_Z[4][oc3] = float(row[0])
                
                if i in [0, 1, 2, 3]:
                # these include extra colums X and Y, filters from 5
    
                    for k in range(5, 15):
                        data_Z[k][oc3] = float(row[k+1])
    
                elif i in [4, 5, 6, 7]:
                # these exclude extra colums X and Y, filters from 3
     
                    for k in range(5, 15):
                        data_Z[k][oc3] = float(row[k-1])  
                        
                else:
                    print("lemon you melon")
                
                oc3 += 1
                
### ANOTHER BASIC CHECK TO CONFIRM NO MAJOR ERRORS ###

# this is checking the ID value from Z and A catalogs match per row
for i in range(len(data_Z[0])): 
    if data_Z[1][i] != data_Z[4][i]:
        print("You're incompetent, Lester")

    
### REMOVING KNOWN SPURIOUS ROWS IN INPUT FILES, redshift of -1 or 99 ###

delete_list = []

for i in range(object_count_z):
    if data_Z[2][i] < 0 or data_Z[2][i] == 99:
        delete_list.append(int(data_Z[3][i]))
        
data_Z = np.delete(data_Z, delete_list, axis=1)
#data_Z = np.delete(data_Z, 3, axis=0)


### PLOTTING number per redshift ###

#plt.yscale('log')
plt.hist(data_Z[2], bins=np.linspace(0, 10.00000001, 11)) # this puts 10.0 max value in bin 9-10.
plt.show()


### PLOTTING number per magnitude ###

plt.hist(data_Z[13], bins=np.linspace(15, 33, 40))
plt.show()


### PER CATALOG Same thing as below using far simpler code, plotting spitzer 3.6 per catalog
# MUCH NICER 

for i in range(8):
    
    # this is data_Z, column 13 (spitzer 3.6), per i (catalog)
    data_Z_13_i = data_Z[13][data_Z[0] == i]
    plt.hist(data_Z_13_i, bins=np.linspace(15, 33, 40))
plt.show()



### PER REDSHIFT Plotting spitzer 3.6 (13) per redshift bin (combined catalog)

for i in range(10):
    
    plt.hist(data_Z[13][ (data_Z[2] > i) & (data_Z[2] <= i+1) ], bins=np.linspace(15, 33, 40))  
    #plt.hist(data_Z[13][ np.trunc(data_Z[2]) == i ], bins=np.linspace(15, 33, 40))  
    print(i, i+1)
    print(len(data_Z[13][ np.trunc(data_Z[2]) == i ]))
plt.show()
    
    


### Same for H160 (11), number per magnitude, plotted at each redshift. 

for i in range(10):
    
    plt.hist(data_Z[11][ np.trunc(data_Z[2]) == i ], bins=np.linspace(15, 33, 40))  
    print(i, i+1)
    print(len(data_Z[11][ np.trunc(data_Z[2]) == i ]))
plt.show()














### Backwards approach to get individual data plotted given strange distribution
'''
for i in range(8):
    print(np.count_nonzero(data_Z[0] == i))
print(len(data_Z[0]))

start	length
0	      3527
3527	  3411
6938	  3388
10326	  3732
14058	  4068
18126	  3448
21574	  4351
25925	  3403
29328	
'''
'''
data_Z_13_0 = data_Z[13][0:3527]
data_Z_13_1 = data_Z[13][3527:6938]
data_Z_13_2 = data_Z[13][6938:10326]
data_Z_13_3 = data_Z[13][10326:14058]
data_Z_13_4 = data_Z[13][14058:18126]
data_Z_13_5 = data_Z[13][18126:21574]
data_Z_13_6 = data_Z[13][21574:25925]
data_Z_13_7 = data_Z[13][25925:29328]

print(len(data_Z_13_0))
print(len(data_Z_13_1))
print(len(data_Z_13_2))
print(len(data_Z_13_3))
print(len(data_Z_13_4))
print(len(data_Z_13_5))
print(len(data_Z_13_6))
print(len(data_Z_13_7))

plt.hist(data_Z_13_0, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_1, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_2, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_3, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_4, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_5, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_6, bins=np.linspace(15, 33, 40))
plt.show()

plt.hist(data_Z_13_7, bins=np.linspace(15, 33, 40))
plt.show()
'''






















