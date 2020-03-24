#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:12:33 2020

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# load the templates
# =============================================================================
load = 100000

alpha = np.load('alpha.npy')[:load]
beta = np.load('beta.npy')[:load]
tau = np.load('tau.npy')[:load]
mass = np.load('mass.npy')[:load]
A = np.load('A.npy')[:load]
sfr = np.load('sfr.npy')[:load]

subset = sfr >= 0

alpha = alpha[subset]
beta = beta[subset]
tau = tau[subset]
mass = mass[subset]
A = A[subset]
sfr = sfr[subset]

print(len(sfr))

# =============================================================================
# sample mass
# =============================================================================
nObj = 100

minMass = 7.
maxMass = 12.

#assign masses drawn from broad Guassian distribution
Marr = np.random.normal(size=nObj,scale=0.5,loc=(minMass+maxMass)/2)

#rechoose a mass if not within minMass and maxMass
while np.max(Marr) > maxMass or np.min(Marr) < minMass:
    tempIdx = np.where((Marr > maxMass) | (Marr < minMass))[0]
    print(tempIdx)
    Marr[tempIdx] = np.random.normal(size=len(tempIdx),scale=0.5,loc=(minMass+maxMass)/2)
    
# =============================================================================
# calculate sfr from linear relation + scatter
# =============================================================================

intercept = -8.
slope = 1.
scatter = 0.3

Sfr = slope*Marr+intercept+np.random.normal(size=nObj,scale=scatter)


# =============================================================================
# PLOT 1
# =============================================================================
cb = tau

xlin = np.array([7, 12])
ylin = slope*xlin + intercept

mock_fit = np.polyfit(Marr, Sfr, deg=1)

plt.figure(figsize=(12, 10))
plt.xlim(7, 12)
plt.ylim(-1, 3)

plt.scatter(np.log10(mass), np.log10(sfr), s=5, c=cb) # templates
plt.plot(xlin, ylin, color='k') # straight line - input relation
plt.colorbar()
plt.scatter(Marr, Sfr) # mock sample
plt.plot(xlin, mock_fit[0]*xlin + mock_fit[1], color='r') # straight line fit - mock sample

plt.show()

# =============================================================================
# find nearest point for single example
# =============================================================================

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

nodes = np.empty([len(mass), 2])

for i in range(len(mass)):
    nodes[i, 0] = np.log10(mass[i])
    nodes[i, 1] = np.log10(sfr[i])

node = [9.8, 0.4]

ind = closest_node(node, nodes)


# =============================================================================
# PLOT 2
# =============================================================================
cb = tau

xlin = np.array([7, 12])
ylin = slope*xlin + intercept

plt.figure(figsize=(12, 10))

plt.scatter(np.log10(mass), np.log10(sfr), s=5, c=cb) # templates
plt.scatter(np.log10(mass)[ind], np.log10(sfr)[ind], s=100, zorder=2, color='r') # nearest template
plt.plot(xlin, ylin, color='k')
plt.plot(node[0], node[1], marker='x', mew=50, zorder=1, color='k', ms=5) # single chosen point


plt.xlim(7, 12)
plt.ylim(-1, 3)

plt.show()











