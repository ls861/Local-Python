#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 18:15:34 2020

@author: lester
"""



import numpy as np
import matplotlib.pyplot as plt

lm = np.load(sbf+'lm_chain_{}x_{}_{}_z_dependence.npy'.format(str(nChains), str(minIter), comment2))



print(lm['xi'])

print(lm['eta'])

plt.hist(lm['xi'][:,0])
plt.hist(lm['eta'][:,0])