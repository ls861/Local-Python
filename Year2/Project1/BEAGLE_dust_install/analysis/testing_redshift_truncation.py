#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 22:33:08 2021

@author: lester

https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/

"""


import numpy as np
from scipy.stats import norm, uniform
import matplotlib.pyplot as plt


# target distribution


# z1 = 0.0
# z2 = np.inf
# x_curr = 3.5
# sigma = 0.5
# def pi_x(x):
#     s = x
#     return (s) * np.exp(-(s)) 

# z1 = 3.0
# z2 = np.inf
# x_curr = 3.5
# sigma = 0.5
# def pi_x(x):
#     s = x - 3.0
#     return (s) * np.exp(-(s)) # normalisation

z1 = 3.0
z2 = 4.0
x_curr = 3.5
sigma = 0.5
def pi_x(x):
    s = 10.0 * (x - 3.0)
    return (s) * np.exp(-(s)) / 0.0999501 # normalisation
    

x_arr = []
x = np.linspace(z1, z2, 1000)






for i in range(3000):
    if i % 10000 == 0:
        print(i)
    x_new = norm.rvs(loc=x_curr, scale=sigma)
    while x_new < z1 or x_new > z2:
        x_new = norm.rvs(loc=x_curr, scale=sigma)

    A = pi_x(x_new) / pi_x(x_curr)
    # A = A * (norm.cdf(x_curr, loc=0, scale=sigma) / norm.cdf(x_new, loc=0, scale=sigma))
    # A = A * (norm.cdf(x_curr-z1, loc=0, scale=sigma) / norm.cdf(x_new-z1, loc=0, scale=sigma))
    A = A * ((norm.cdf(x_curr-z1, loc=0, scale=sigma)-norm.cdf(x_curr-z2, loc=0, scale=sigma)) / (norm.cdf(x_new-z1, loc=0, scale=sigma)-norm.cdf(x_new-z2, loc=0, scale=sigma)))
    accProb = min(1, A)
    
    x_curr_old = x_curr
    
    u = uniform.rvs()
    if u <= accProb:
        x_curr = x_new

    x_arr.append(x_curr)

plt.hist(x_arr, bins=70, density=True)
plt.plot(x, pi_x(x))
plt.show()

# plt.plot(x_arr)
# plt.show()

# print((norm.cdf(x_curr-z1)-norm.cdf(x_curr-z2)))
# print((norm.cdf(x_new-z1)-norm.cdf(x_new-z2)))
# print(x_curr_old, x_new)
# print(((norm.cdf(x_curr_old-z1)-norm.cdf(x_curr_old-z2)) / (norm.cdf(x_new-z1)-norm.cdf(x_new-z2))))

#%%
import numpy as np
from scipy.stats import norm, uniform
import matplotlib.pyplot as plt

z1 = 0.0
z2 = np.inf
x_curr = 3.5
sigma = 0.5
def pi_x2(x):
    s = x
    return (s) * np.exp(-(s)) 

x_arr = []
x = np.linspace(0, 10, 1000)

# let proposal be a normal distribution
x_curr = 5.0
sigma = 1.0

for i in range(3000):
    x_new = norm.rvs(loc=x_curr, scale=sigma)
    while x_new < 0.0:
        x_new = norm.rvs(loc=x_curr, scale=sigma)

    A = pi_x2(x_new) / pi_x2(x_curr)
    A = A * (norm.cdf(x_curr) / norm.cdf(x_new))
    accProb = min(1, A)
    
    u = uniform.rvs()
    if u <= accProb:
        x_curr = x_new

    x_arr.append(x_curr)

plt.hist(x_arr, bins=70, density=True)
plt.plot(x, pi_x2(x))
plt.show()

plt.plot(x_arr)
plt.show()

#%%
# =============================================================================
# deciding whether I can have prior on top of truncation...
# https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5B10.0+*+%5C%2840%29x+-+3.0%5C%2841%29+*+exp%5C%2840%29-10.0+*+%5C%2840%29x+-+3.0%5C%2841%29%5C%2841%29%2C%7Bx%2C3.1%2C4.0%7D%5D
# =============================================================================


# target distribution

z1 = 3.0
z2 = 4.0
x_curr = 3.5
sigma = 0.5
def pi_x(x):
    s = 10.0 * (x - 3.0)
    return (s) * np.exp(-(s)) / 0.0735259 # normalisation from x=3.1 - 4.0
    

x_arr = []
x = np.linspace(z1, z2, 1000)


for i in range(300000):
    if i % 1000 == 0:
        print(i)
    x_new = norm.rvs(loc=x_curr, scale=sigma)
    while x_new < z1 or x_new > z2:
        x_new = norm.rvs(loc=x_curr, scale=sigma)

    A = pi_x(x_new) / pi_x(x_curr)
    # A = A * (norm.cdf(x_curr, loc=0, scale=sigma) / norm.cdf(x_new, loc=0, scale=sigma))
    # A = A * (norm.cdf(x_curr-z1, loc=0, scale=sigma) / norm.cdf(x_new-z1, loc=0, scale=sigma))
    A = A * ((norm.cdf(x_curr-z1, loc=0, scale=sigma)-norm.cdf(x_curr-z2, loc=0, scale=sigma)) / (norm.cdf(x_new-z1, loc=0, scale=sigma)-norm.cdf(x_new-z2, loc=0, scale=sigma)))
    accProb = min(1, A)
    
    x_curr_old = x_curr
    
    u = uniform.rvs()
    if u <= accProb:
        
        # PRIOR ON TOP 
        if x_new > 3.1:
            x_curr = x_new

    x_arr.append(x_curr)

plt.hist(x_arr, bins=70, density=True)
plt.plot(x, pi_x(x))
plt.show()

# plt.plot(x_arr)
# plt.show()

# print((norm.cdf(x_curr-z1)-norm.cdf(x_curr-z2)))
# print((norm.cdf(x_new-z1)-norm.cdf(x_new-z2)))
# print(x_curr_old, x_new)
# print(((norm.cdf(x_curr_old-z1)-norm.cdf(x_curr_old-z2)) / (norm.cdf(x_new-z1)-norm.cdf(x_new-z2))))







