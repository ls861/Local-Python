#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:40:15 2020

@author: lester
"""

import numpy as np

x = np.linspace(1000000, 11E8, 7)

print(len(x))
print(x)

print(13*8)

print(10%4)


print(np.log10(x))



for i in range(7):
    print(np.log10(x[i]))
