#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:00:57 2020

@author: lester
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('figure', figsize=(4,4))

def best_fit(x, poly):
    
    y = 0
    for i in range(len(poly)):
        y += (x**(len(poly)-(1+i))) * poly[i]
#    y = poly[0]*x + poly[1]
#    y = poly[0]*(x**4) + poly[0]*(x**3) + poly[0]*(x**2) + poly[0]*(x**1) + poly[1]
    return y

# =============================================================================
# CREATING A MODEL
# =============================================================================
    
run = '004'

runtime = np.load('DE_108_{}_runtime.npy'.format(run))
log_H160_a = np.load('DE_108_{}_log_H160_a.npy'.format(run))


poly = np.polyfit(log_H160_a, runtime, 2)

model_runtime = best_fit(log_H160_a, poly)

print(max(log_H160_a))
print(min(runtime))
print(sum(runtime))

print(min(model_runtime))
print(sum(model_runtime))


plt.scatter(log_H160_a, runtime)
plt.scatter(log_H160_a, model_runtime, color='r')
plt.show()

# =============================================================================
# APPLYING MODEL TO ALL OF ASTRODEEP
# =============================================================================

D = np.load('astrodeep_rawfile.npy')
print(len(D))


print(D['b_H160'][D['b_H160'] > 500])

print(min  (D['b_H160'][D['b_H160'] > 0]    )    )
# note, there are -ve values, 0 values, and lowest other value is 7.881e-06
# adding 1e-6 as usual to prevent zeros



D_log_H160_a = np.log10(abs(D['b_H160'])+1e-6)

                   
plt.scatter(D_log_H160_a, best_fit(D_log_H160_a, poly), color='g')            
plt.scatter(log_H160_a, runtime)
plt.scatter(log_H160_a, model_runtime, color='r')
plt.show()               

plt.hist(best_fit(D_log_H160_a, poly))
plt.show()

print(sum(best_fit(D_log_H160_a, poly)))
print(min(best_fit(D_log_H160_a, poly)))
print(max(best_fit(D_log_H160_a, poly)))

print('RUN: ' + str(run))
print('H160 min (uJy)', '# objects', 'runtime (h)', 'inc max runtime')
for i in [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]:


    # enforcing a hard upper limit
    upper = max(best_fit(log_H160_a, poly)) # 3.1643
    upper_arr = np.full(len(D_log_H160_a), upper)
    
    print(10**i, len(best_fit(D_log_H160_a, poly)[D_log_H160_a > i]), sum(best_fit(D_log_H160_a, poly)[D_log_H160_a > i]), sum(np.minimum(best_fit(D_log_H160_a, poly), upper_arr)[D_log_H160_a > i]))
        
# =============================================================================
# Characterising ASTRODEEP
# =============================================================================
    
       
minus = D['b_H160'][D['b_H160'] < 0]
zero = D['b_H160'][D['b_H160'] == 0]
plus = D['b_H160'][D['b_H160'] > 0]

print(len(minus), len(zero), len(plus))




plt.hist(minus, bins=100)
plt.yscale('log')
plt.show()

plt.hist(plus, bins=100)
plt.yscale('log')
plt.show()  
    
plt.hist(plus, bins=100, range=(0, 2000))
plt.yscale('log')
plt.show()  
    
    
test = np.sort(D, order='b_H160')[:-10:-1]
print(test['b_H160'])
print(test['field'])
print(test['ID'])   

test = np.sort(D, order='b_H160')[:10]
print(test['b_H160'])
print(test['field'])
print(test['ID'])   
    
    

mpl.rcdefaults()
























