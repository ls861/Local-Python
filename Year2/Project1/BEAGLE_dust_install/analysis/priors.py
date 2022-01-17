#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 12:42:57 2021

@author: lester
"""

import numpy as np
import matplotlib.pyplot as plt

n = 1000000
mass = 10**np.random.uniform(5, 12, n) # 5 12
msa  = 10**np.random.uniform(6, 10, n) # 6 10
tau = 10**np.random.uniform(7, 10.5, n) # 7 10.5
tau_i = 1/np.random.uniform(1.0/(10**10), 1.0/(10**7.0), n)



log10_mass = np.log10(mass)
log10_msa = np.log10(msa)
log10_tau = np.log10(tau)
log10_tau_i = np.log10(tau_i)
'''
plt.hist(log10_mass, histtype='step')
plt.hist(log10_msa, histtype='step')
plt.hist(log10_tau, histtype='step')
plt.show()

plt.hist(log10_tau, histtype='step', density='True', label='uniform on log tau', color='r')
plt.hist(log10_tau_i, histtype='step', density='True', label='uniform on 1/tau', color='b')
plt.legend()
plt.show()
'''
x_tmp = np.array((10**5, 10**12))
log10_x_tmp = np.log10(x_tmp)



log10_sfr_c = log10_mass - log10_msa

norm_denom = (-tau*np.exp(-msa/tau)*(tau+msa)) + np.power(tau,2)
log10_norm = log10_mass - np.log10(norm_denom)
exp_term = -msa/tau
exp_idx = (exp_term <= -100.0)
exp_term[exp_idx] = -100.0
temp_sfr = log10_norm + log10_msa + np.log10(np.exp(exp_term))
log10_sfr_de = np.where(temp_sfr < -30.0, -30.0, temp_sfr)
idx = log10_msa < log10_tau

norm_denom_i = (-tau_i*np.exp(-msa/tau_i)*(tau_i+msa)) + np.power(tau_i,2)
log10_norm_i = log10_mass - np.log10(norm_denom_i)
exp_term_i = -msa/tau_i
exp_idx_i = (exp_term_i <= -100.0)
exp_term_i[exp_idx_i] = -100.0
temp_sfr_i = log10_norm_i + log10_msa + np.log10(np.exp(exp_term_i))
log10_sfr_de_i = np.where(temp_sfr_i < -30.0, -30.0, temp_sfr_i)
idx_i = log10_msa < log10_tau_i


# constant
idx_tmp = abs(log10_sfr_c) < 5
plt.hist2d(log10_mass[idx_tmp], log10_sfr_c[idx_tmp], bins=(50, 50))
plt.plot(log10_x_tmp, log10_x_tmp - 6, c='r', label='constant 10p6')
plt.plot(log10_x_tmp, log10_x_tmp - 6 + np.log10(2), c='r', linestyle='dotted', label='linear 10p6') # t = 6
plt.plot(log10_x_tmp, log10_x_tmp - 10, c='r', label='constant 10p10') # t = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 - np.log10(np.exp(1.0)-2.0), c='r', linestyle='dashed', label='ridge') # ridge if t = tau = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 + np.log10(2), c='r', linestyle='dotted', label='linear 10p10') # t = 10
plt.ylim(-5, 5)
plt.show()

# DE uniform log tau
idx_tmp = abs(log10_sfr_de) < 5
plt.hist2d(log10_mass[idx_tmp], log10_sfr_de[idx_tmp], bins=(50, 50))
plt.plot(log10_x_tmp, log10_x_tmp - 6, c='r', label='constant 10p6')
plt.plot(log10_x_tmp, log10_x_tmp - 6 + np.log10(2), c='r', linestyle='dotted', label='linear 10p6') # t = 6
plt.plot(log10_x_tmp, log10_x_tmp - 10, c='r', label='constant 10p10') # t = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 - np.log10(np.exp(1.0)-2.0), c='r', linestyle='dashed', label='ridge') # ridge if t = tau = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 + np.log10(2), c='r', linestyle='dotted', label='linear 10p10') # t = 10
plt.ylim(-5, 5)
plt.show()


# DE uniform 1/tau
idx_tmp = abs(log10_sfr_de_i) < 5
plt.hist2d(log10_mass[idx_tmp], log10_sfr_de_i[idx_tmp], bins=(50, 50))
plt.plot(log10_x_tmp, log10_x_tmp - 6, c='r', label='constant 10p6')
plt.plot(log10_x_tmp, log10_x_tmp - 6 + np.log10(2), c='r', linestyle='dotted', label='linear 10p6') # t = 6
plt.plot(log10_x_tmp, log10_x_tmp - 10, c='r', label='constant 10p10') # t = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 - np.log10(np.exp(1.0)-2.0), c='r', linestyle='dashed', label='ridge') # ridge if t = tau = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 + np.log10(2), c='r', linestyle='dotted', label='linear 10p10') # t = 10
plt.ylim(-5, 5)
plt.show()


# DE uniform log tau RISING ONLY
idx_tmp = abs(log10_sfr_de[idx]) < 5
plt.hist2d(log10_mass[idx][idx_tmp], log10_sfr_de[idx][idx_tmp], bins=(50, 50))
plt.plot(log10_x_tmp, log10_x_tmp - 6, c='r', label='constant 10p6')
plt.plot(log10_x_tmp, log10_x_tmp - 6 + np.log10(2), c='r', linestyle='dotted', label='linear 10p6') # t = 6
plt.plot(log10_x_tmp, log10_x_tmp - 10, c='r', label='constant 10p10') # t = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 - np.log10(np.exp(1.0)-2.0), c='r', linestyle='dashed', label='ridge') # ridge if t = tau = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 + np.log10(2), c='r', linestyle='dotted', label='linear 10p10') # t = 10
plt.ylim(-5, 5)
plt.show()


# DE uniform log tau RISING ONLY
idx_tmp = abs(log10_sfr_de[idx]) < 5
plt.scatter(log10_mass[idx][idx_tmp], log10_sfr_de[idx][idx_tmp], s=0.01)
plt.plot(log10_x_tmp, log10_x_tmp - 6, c='r', label='constant 10p6')
plt.plot(log10_x_tmp, log10_x_tmp - 6 + np.log10(2), c='r', linestyle='dotted', label='linear 10p6') # t = 6
plt.plot(log10_x_tmp, log10_x_tmp - 10, c='r', label='constant 10p10') # t = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 - np.log10(np.exp(1.0)-2.0), c='r', linestyle='dashed', label='ridge') # ridge if t = tau = 10
plt.plot(log10_x_tmp, log10_x_tmp - 10 + np.log10(2), c='r', linestyle='dotted', label='linear 10p10') # t = 10
plt.ylim(-5, 5)
plt.show()



#%%
# fixed mass

n = 1000000
mass = 10**np.random.uniform(6, 6, n) # 5 12
msa  = 10**np.random.uniform(6, 10, n) # 6 10
tau = 10**np.random.uniform(7, 10, n) # 7 10.5
tau_i = 1/np.random.uniform(1.0/(10**10), 1.0/(10**7.0), n)

log10_mass = np.log10(mass)
log10_msa = np.log10(msa)
log10_tau = np.log10(tau)
log10_tau_i = np.log10(tau_i)

x_tmp = np.array((10**5, 10**12))
log10_x_tmp = np.log10(x_tmp)

log10_sfr_c = log10_mass - log10_msa

norm_denom = (-tau*np.exp(-msa/tau)*(tau+msa)) + np.power(tau,2)
log10_norm = log10_mass - np.log10(norm_denom)
exp_term = -msa/tau
exp_idx = (exp_term <= -100.0)
exp_term[exp_idx] = -100.0
temp_sfr = log10_norm + log10_msa + np.log10(np.exp(exp_term))
log10_sfr_de = np.where(temp_sfr < -30.0, -30.0, temp_sfr)
idx = log10_msa < log10_tau

norm_denom_i = (-tau_i*np.exp(-msa/tau_i)*(tau_i+msa)) + np.power(tau_i,2)
log10_norm_i = log10_mass - np.log10(norm_denom_i)
exp_term_i = -msa/tau_i
exp_idx_i = (exp_term_i <= -100.0)
exp_term_i[exp_idx_i] = -100.0
temp_sfr_i = log10_norm_i + log10_msa + np.log10(np.exp(exp_term_i))
log10_sfr_de_i = np.where(temp_sfr_i < -30.0, -30.0, temp_sfr_i)
idx_i = log10_msa < log10_tau_i


plt.hist(log10_sfr_de[log10_sfr_de>-5], bins=20, density=True, histtype='step', color='r', label='uniform log tau')
plt.hist(log10_sfr_de_i[log10_sfr_de_i>-5], bins=20, density=True, histtype='step', color='b', label='uniform 1/tau')
plt.legend(loc='upper left')
plt.show()



#%%



from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
z_med_spe = np.linspace(0.0, 8.0, 100)
t_med_spe = cosmo.age(z_med_spe).value

t_now = cosmo.age(0.0).value
print(t_now)

plt.plot(z_med_spe, t_med_spe, label='age of universe')
plt.plot(z_med_spe, t_now-t_med_spe, label='lookback time')
plt.xlabel('redshift')
plt.ylabel('age of universe / Gyr')
plt.legend()
plt.show()
















