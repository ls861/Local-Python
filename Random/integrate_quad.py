#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:34:05 2020

@author: lester
"""

from scipy.integrate import quad
import numpy as np






A      = 5
t0     = 6
tau    = 1
t      = 10

#t      = np.arange(t0+0.01, 15, 0.1)

sfr             = A * np.heaviside(t-t0, 1) * (1/t) * np.exp( -(((np.log(t-t0))**2) / (2*(tau**2))) )
integrand       = lambda T: np.heaviside(T-t0, 1) * (1/T) * np.exp( -(((np.log(T-t0))**2) / (2*(tau**2))) )






i = quad(integrand, 6, 8)
print(i)

m = 10**7

A = m / i[0]
print(A)







