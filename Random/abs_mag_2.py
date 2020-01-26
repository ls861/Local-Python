#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:44:00 2020

@author: lester
"""

import numpy as np

m = 23.1296
z = 7.1385
Dl = 70573.1 * (10**6) #pc

bracket = (1/(1 + z)) * ( (Dl / 10) ** 2 )

M = m - (2.5* np.log10(bracket))

print(M)












