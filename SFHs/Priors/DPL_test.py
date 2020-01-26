#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:05:06 2020

@author: lester
"""

import numpy as np


A = 2.
alpha = 10.

tau = 1E9
t = 5E8
k = 100.



beta_min = -(np.log((2*k) - ((t/tau)**(alpha)))) / np.log(t/tau)




beta = beta_min + 0.1


sfr = A / (((t/tau)**alpha)+((t/tau)**-beta))
sfr_tau = A / (((tau/tau)**alpha)+((tau/tau)**-beta))


print(sfr, sfr_tau)
print(beta_min)

#('NANANANAN', nan, 1157523536.5616436, 1013011483.8305516, 83.18827835263626)
#('NANANANAN', nan, 1204936679.3669977, 1063061402.2575912, 90.63399003839977)


t_arr = 1157523536.5616436
tau_arr = 1013011483.8305516

t_arr = 60
tau_arr = 20.
alpha_arr = 5.

beta_max    = -(np.log10((2*k) - ((t_arr/tau_arr)**(alpha_arr)))) / np.log10(t_arr/tau_arr)

print(beta_max)



# beta = -inf when beta = tau.


