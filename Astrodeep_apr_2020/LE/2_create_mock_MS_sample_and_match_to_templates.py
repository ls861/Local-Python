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
load = 1000000 # max 1,000,000

msa = np.load('msa.npy')[:load]
tau = np.load('tau.npy')[:load]
tau_exp = np.load('tau_exp.npy')[:load]
mass = np.load('mass.npy')[:load]
A = np.load('A.npy')[:load]
sfr = np.load('sfr.npy')[:load]

subset = sfr >= 0

msa = msa[subset]
tau = tau[subset]
tau_exp = tau_exp[subset]
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

## ORIGINAL
#intercept = -8.
#slope = 1.
#scatter = 0.3

ageUniv = 3228839870.9122815 # from previous python script, z=2, {'o_M_0':0.3,'o_lambda_0':0.7,'h':0.7}
slope = 0.84 - 0.026*ageUniv*1e-9 # 0.7560501633562806
intercept = -(6.51 - 0.11*ageUniv*1e-9) # -6.1548276141996485
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

# =============================================================================
# find nearest point for entire MS sample
# =============================================================================

closest_ind = np.empty(nObj, dtype=int)

for i in range(nObj):
    closest_ind[i] = closest_node([Marr[i], Sfr[i]], nodes)

plt.figure(figsize=(12, 10))
plt.scatter(Marr, Sfr, marker='x', zorder=1)
plt.scatter(np.log10(mass)[closest_ind], np.log10(sfr)[closest_ind], zorder=0)





# =============================================================================
#np.save('Marr', Marr)
#np.save('sfr_mock', Sfr)
#np.save('closest_ind', closest_ind)
# =============================================================================




















