#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 00:45:00 2020

@author: lester
"""




import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm

import matplotlib as mpl
#mpl.rcParams.update(mpl.rcParamsDefault)


# =============================================================================
# investigating multiple slices and then multple k values
# =============================================================================


count = 100000

xlow = 1.0
xhigh = 20.0

x = np.linspace(xlow, xhigh, count)
x = np.random.normal(9.3, 0.4, count)

alpha = -6.2
beta = 0.8
sig_0 = 0.3

xi_min = 8.5
xi_max = 10.0
k = 2 # factor higher sigma at min mass compared to max mass

sigma = sig_0 * ( ((1.0-k)*(x-xi_max) / (xi_max - xi_min)) + 1)
#print(sigma)

y4 = alpha + (beta*x) + np.random.normal(0, abs(sigma), len(x))

plt.hist2d(x, y4, bins=50)
plt.show()



xlow = 8.0
xhigh = 10.5
ylow = -0.5
yhigh = 3
siglow = sig_0 * ( ((1.0-k)*(xlow-xi_max) / (xi_max - xi_min)) + 1)
sighigh = sig_0 * ( ((1.0-k)*(xhigh-xi_max) / (xi_max - xi_min)) + 1)

plt.figure(figsize=(10, 6))
plt.title('Main Sequence', size = 20)
plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)
plt.scatter(x, y4, marker='x')
plt.plot((xlow, xhigh), (alpha+xlow*beta, alpha+xhigh*beta), color='k')
plt.plot((xi_min, xi_min), (ylow, yhigh), color='r')
plt.plot((xi_max, xi_max), (ylow, yhigh), color='r')
plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='r')
plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='r')
plt.show()

plt.figure(figsize=(10, 6))
plt.title('Main Sequence', size = 20)
plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
#plt.xlim(7.5, 11)
#plt.ylim(-4, 4)
plt.hist2d(x, y4, bins=50, range=((xlow, xhigh),(ylow, yhigh)), norm=mcolors.LogNorm(), cmap=plt.cm.viridis)
#plt.colorbar()
cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax = plt.axes()
ax.set_facecolor(rgba)
plt.plot((xlow, xhigh), (alpha+xlow*beta, alpha+xhigh*beta), color='k')
#plt.plot((xi_min, xi_min), (ylow, yhigh), color='r')
#plt.plot((xi_max, xi_max), (ylow, yhigh), color='r')
plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='r')
plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='r')
plt.show()


# for poster
plt.figure(figsize=(10, 6))
lw = 4
#plt.title('Main Sequence', size = 20)
plt.xlabel(r'$\mathrm{log}(m_{tot}/M_{\odot})$', size = 16)
plt.ylabel(r'$\mathrm{log}(\Psi / M_{\odot} yr^{-1})$', size = 16)
plt.hist2d(x, y4, bins=50, range=((xlow, xhigh),(ylow, yhigh)), cmap=plt.cm.viridis)
plt.colorbar(label='Count')
cmap = cm.get_cmap('viridis')
rgba = cmap(0)
ax = plt.axes()
ax.set_facecolor(rgba)
plt.plot((xlow, xhigh), (alpha+xlow*beta, alpha+xhigh*beta), color='k', linewidth=lw, label=r'Linear Relation, alpha = -6.2, beta = 0.8')
plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='k', linewidth=lw, linestyle='--', label=r'Intrinsic Scatter, sig0 = 0.3, k = 2')
plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='k', linewidth=lw, linestyle='--')
plt.legend(loc='lower right')
plt.show()


# =============================================================================
# first year report
# =============================================================================






count = 100000

xlow = 1.0
xhigh = 20.0

x = np.linspace(xlow, xhigh, count)
#x = np.random.normal(9.3, 0.4, count)

alpha = -6.2
beta = 0.8
sig_0 = 0.2

xi_min = 8.5
xi_max = 10.0
k = 3 # factor higher sigma at min mass compared to max mass

sigma = sig_0 * ( ((1.0-k)*(x-xi_max) / (xi_max - xi_min)) + 1)
#print(sigma)

y4 = alpha + (beta*x) + np.random.normal(0, abs(sigma), len(x))

plt.hist2d(x, y4, bins=50)
plt.show()



xlow = 8.0
xhigh = 10.5
ylow = -0.5
yhigh = 3
siglow = sig_0 * ( ((1.0-k)*(xlow-xi_max) / (xi_max - xi_min)) + 1)
sighigh = sig_0 * ( ((1.0-k)*(xhigh-xi_max) / (xi_max - xi_min)) + 1)


import matplotlib
matplotlib.rcParams.update({'font.size': 18})
        
plt.figure(figsize=(10, 6))
#plt.title('Main Sequence', size = 20)
plt.xlabel(r'$M_\star$')
plt.ylabel(r'$\Psi$')
plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)
plt.scatter(x, y4, marker='x')
plt.plot((xlow, xhigh), (alpha+xlow*beta, alpha+xhigh*beta), color='k')
plt.plot((xi_min, xi_min), (ylow, yhigh), color='r')
plt.plot((xi_max, xi_max), (ylow, yhigh), color='r')
#plt.plot((xlow, xhigh), (alpha+beta*xlow+siglow, alpha+beta*xhigh+sighigh), color='r')
#plt.plot((xlow, xhigh), (alpha+beta*xlow-siglow, alpha+beta*xhigh-sighigh), color='r')
plt.tight_layout()
#plt.savefig('/Users/lester/Dropbox/PhD/20_Summer/First Year Report/RawFigs/320_mass_dep_model.png')
plt.show()














