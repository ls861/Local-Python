#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:30:23 2020

@author: lester
"""


import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, np.pi, 10000)

x = np.cos(14*t)
y = np.sin(26*t)

plt.plot(x, y)






print(np.pi)