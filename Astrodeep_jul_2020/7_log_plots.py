#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 11:48:38 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt




x = np.linspace(0, 100, 1000)
y = x

plt.plot(x, y)
plt.show()

plt.plot(np.log10(x), y)
plt.show()

plt.plot(x, np.log10(y))
plt.show()

plt.plot(np.log10(x), np.log10(y))
plt.show()



x = np.linspace(0, 100, 1000)
y = x**2

plt.plot(x, y)
plt.show()

plt.plot(np.log10(x), y)
plt.show()

plt.plot(x, np.log10(y))
plt.show()

plt.plot(np.log10(x), np.log10(y))
plt.show()



x = np.linspace(0, 100, 1000)
y = x**2 +10

plt.plot(x, y)
plt.show()

plt.plot(np.log10(x), y)
plt.show()

plt.plot(x, np.log10(y))
plt.show()

plt.plot(np.log10(x), np.log10(y))
plt.show()
























