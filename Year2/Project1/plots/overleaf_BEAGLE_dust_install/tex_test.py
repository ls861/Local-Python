#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 21:52:59 2021

@author: lester
"""



import matplotlib
matplotlib.__version__


import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# import matplotlib.pyplot as plt
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})

# # for Palatino and other serif fonts use:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Palatino"],
# })

plt.plot((0,0), (0,0))
plt.ylabel(r'$\mu$')
plt.show()




