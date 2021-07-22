#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:27:23 2020

@author: lester
"""



import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0.1, 10.0, 100)


y = np.array([x, 3*x, x+10, (3*x)+10, x**2, x**3, 3.0*(x**4.0)])

y = np.array([x, x**2, x**3, 10*x, 10*(x**2), 10*(x**2)+100])

for i in range(len(y)):
    plt.plot(x, y[i])
plt.show()

for i in range(len(y)):
    plt.plot(x, y[i])
plt.xscale('log')
plt.show()

for i in range(len(y)):
    plt.plot(x, y[i])
plt.yscale('log')
plt.show()

for i in range(len(y)):
    plt.plot(x, y[i])
plt.xscale('log')
plt.yscale('log')
plt.show()