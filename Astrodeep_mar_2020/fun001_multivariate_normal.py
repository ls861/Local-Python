#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 09:52:43 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt


s=1

mean = (0, 0)
cov = [[0.15, -0.01], [-0.01, 0.15]]
x = np.random.multivariate_normal(mean, cov, 100000)


plt.hist2d(x[:,0], x[:,1], bins=50)
plt.xlim(-s, s)
plt.ylim(-s, s)


















