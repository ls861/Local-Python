#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:05:31 2020

@author: lester
"""

from scipy import stats
import matplotlib.pyplot as plt

mu = 54
sig = 13
k = [1, 5, 10, 20]

size = 10000

Q_arr = []

for k in k:
    Q_arr = []
    for i in range(size):
        dist = stats.norm(mu, sig)  # mean = mu, stdev = sig
        r = dist.rvs(k)             # k random draws
        z = (r-mu) / sig
        Q = sum(z**2)
        Q_arr.append(Q)
    
    
    plt.hist(Q_arr, bins=30, histtype=u'step')

plt.ylim(0, 2000)




























