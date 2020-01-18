#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 17:21:53 2019

@author: lester
"""

import sys
import numpy as np

import matplotlib.pyplot as plt

#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=10)

a = [(1, 2, 3), (1, 2, 3), (1, 2, 3)]
b = [(4, 5, 6), (4, 5, 6), (4, 5, 7)]

c = np.concatenate((a, b), axis=0)
print(c)

c = np.concatenate((a, b), axis=1)
print(c)




print(c.dtype.names)

test = [0, 1, 2, 3, 4, 5, 6, 7]

print(test[0::2])
print(test[1::2])