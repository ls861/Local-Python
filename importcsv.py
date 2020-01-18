#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:57:57 2019

@author: lester
"""

import csv
import matplotlib.pyplot as plt

spz_wl_1 = []
spz_data_1 = []
spz_wl_2 = []
spz_data_2 = []

with open('/Users/lester/BEAGLE/BEAGLE-general/filters/080924ch1trans_full.txt') as csvfile:
    fits = csv.reader(csvfile, delimiter='\t')
    for row in fits:
        spz_wl_1.append(float(row[0]))
        spz_data_1.append(float(row[1]))

with open('/Users/lester/BEAGLE/BEAGLE-general/filters/080924ch2trans_full.txt') as csvfile:
    fits = csv.reader(csvfile, delimiter='\t')
    for row in fits:
        spz_wl_2.append(float(row[0]))
        spz_data_2.append(float(row[1]))
        
        
    
x = [1, 2, 3, 4, 5]
y = [1, 4, 9, 16, 25]

plt.plot(x,y)

print(spz_wl_1)
print(spz_data_1)
plt.plot(spz_wl_1, spz_data_1)

print(spz_wl_1[0])
print(type(spz_wl_1[0]))






