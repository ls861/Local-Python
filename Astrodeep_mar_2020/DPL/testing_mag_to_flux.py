#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 23:00:25 2020

@author: lester
"""

import numpy as np

app_mag = np.array([22.55272, 22.508825, 22.46062, 22.176414, 21.74826, 21.638584, 21.628305, 21.638666, 21.759888, 21.88133])

app_mag = np.array([26.7572, 23.936129, 23.75633, 22.122265, 22.55272])


flux = 10**( (23.9 - app_mag) / 2.5 )


print(flux)













