#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 12:22:03 2020

@author: lester
"""



import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


A = 1
t = np.linspace(0, 1000, 10000)
tau = 300
tau_exp = 40


sfr = A * (t*np.heaviside(tau - t, 0)  + tau*np.exp((tau-t)/tau_exp)*np.heaviside(t - tau, 1))
integrand = lambda T: (T*np.heaviside(tau - T, 0)  + tau*np.exp((tau-T)/tau_exp)*np.heaviside(T - tau, 1))

plt.plot(t, sfr)


i = quad(integrand, 0, 400)

print(i)

