#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:41:48 2020

@author: lester
"""

n=1001

import numpy as np
x = 3*np.random.random(n) # 100 random numbers
q25, q50, q75 = np.percentile(x, [25, 50, 75])

print(q25, q50, q75)

print(sum(x[:250]/250))

x = np.sort(x)

print(x[499], x[500], (x[499]+x[500])/2)
print(x[249], x[250], (x[249]+x[250])/2)


y = sum(x)/n
print(y)
