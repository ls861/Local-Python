#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 22:07:59 2020

@author: lester



HUGELY STREAMLINED THE PROCESS BY ACCIDENT, WORTH PORTING INTO REAL FILE AT SOME POINT
"""



import numpy as np
import matplotlib.pyplot as plt
import pickle


chain = pickle.load(open('/Users/lester/Documents/linmix_files/lm_chain_redshift_test_00q.p','r')) # 4 chains

ndraw = len(chain['alpha_a'])/8
#print(ndraw)

psi={}

#keys=['alpha_a','alpha_b','alpha_c','beta_a','beta_b','beta_c','sig_0','k']
keys=['alpha_a']

for key in keys:
    psi[key] = np.vstack((chain[key][1*ndraw:2*ndraw], chain[key][3*ndraw:4*ndraw], chain[key][5*ndraw:6*ndraw], chain[key][7*ndraw:8*ndraw])).T # c['alpha'] is minIter long

#psi['ximean'] = np.vstack([np.sum(c['pi'][ndraw:2*ndraw] * c['mu'][ndraw:2*ndraw], axis=1) for c in chains]).T
#
#psi['xivar'] = np.vstack([np.sum(c['pi'][ndraw:2*ndraw] * (c['tausqr'][ndraw:2*ndraw] + c['mu'][ndraw:2*ndraw]**2), axis=1) for c in chains]).T - psi['ximean']**2




#keys=['alpha_a','alpha_b','alpha_c','beta_a','beta_b','beta_c','sig_0','k','ximean','xivar']
keys=['alpha_a']
for key in keys:
    psi[key] = np.hstack((psi[key][:ndraw/2],psi[key][-ndraw/2:]))

        
#keys=['alpha_a','alpha_b','alpha_c','beta_a','beta_b','beta_c','sig_0','k','ximean','xivar']
keys=['alpha_a']

n = len(psi[keys[0]]) # quarter of initial iterations, half discarded, then split in half
m = len(psi[keys[0]][0]) # 2x number of chains, usually = 8
#psi_bar_dot_j = {}
#psi_bar_dot_dot = {}
#B = {}
#s_j_sqr = {}
W = {}
var_plus = {}
Rhat = {}

test_psi_bar_dot_j = {}
test_psi_bar_dot_dot = {}
test_B = {}
test_s_j_sqr = {}
test_W = {}

for key in keys:
#    psi_bar_dot_j[key] = np.empty(m, dtype=float)
    test_psi_bar_dot_j[key] = np.empty(m, dtype=float)
    test_psi_bar_dot_j[key] = np.mean(psi[key], axis=0) # 8x avg values of each chain
#    for j in range(m):
#        psi_bar_dot_j[key][j] = sum(psi[key][:,j])/n

#    psi_bar_dot_dot[key] = sum(test_psi_bar_dot_j[key])/m
    test_psi_bar_dot_dot[key] = np.mean(test_psi_bar_dot_j[key]) # 1x average of averages

#    B[key] = (n/(m-1.0)) * sum((test_psi_bar_dot_j[key]-test_psi_bar_dot_dot[key])**2) # between sequence
    test_B[key] = (m*n/(m-1.0)) * np.mean((test_psi_bar_dot_j[key]-test_psi_bar_dot_dot[key])**2) # between sequence   
    
#    s_j_sqr[key] = np.empty(m, dtype=float)
#    for j in range(m):
#        s_j_sqr[key][j] = (1.0/(n-1.0)) * sum((psi[key][:,j]-test_psi_bar_dot_j[key][j])**2)
    test_s_j_sqr[key] = np.empty(m, dtype=float)
    test_s_j_sqr[key] = np.sum((psi[key]-test_psi_bar_dot_j[key])**2, axis=0) / (n-1.0)
    
    
#    W[key] = sum(test_s_j_sqr[key])/m # within sequence
    test_W[key] = np.mean(test_s_j_sqr[key]) # within sequence
    
    var_plus[key] = (((n-1.0)/n)*test_W[key]) + (test_B[key]/n)
    
    Rhat[key] = np.sqrt(var_plus[key]/test_W[key])

print(test_B[key]/n, var_plus[key], test_W[key], Rhat[key])

print(var_plus[key]/test_W[key])
print((var_plus[key]/test_W[key])**0.5)
print(1.01**0.5)
print(n)

print((2499./2500.)+(test_B[key]/n)/test_W[key])
print(((2499./2500.)+(test_B[key]/n)/test_W[key])**0.5)

print(2499./2500.)



plt.figure(figsize=(10,10))
for i in range(8):
    plt.plot(psi[key][:,i])
plt.show()







