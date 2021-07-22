#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:30:49 2020

@author: lester
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate




age = 1553026000.2
tau = 14858538000.0
mass = 1e9

time = np.linspace(0, 1e10, 1000)
sfr = time*np.exp(-time/tau)


plt.plot(time, sfr)
plt.plot((age, age), (0, max(sfr)))
plt.show()



sfr_function = lambda t: t*np.exp(-t/tau)
norm_denom_real = integrate.quad(sfr_function, 0, age)[0]






norm_denom = ( (-tau*np.exp(-age/tau)*(tau+age) )+np.power(tau,2))

norm_denom_emma = tau*np.exp(-age/tau)*(tau+age)   +np.power(tau,2)

norm_denom_new = -( (tau*np.exp(-age/tau)*(tau+age) )+np.power(tau,2))

norm_denom_wolf = tau * (tau - (np.exp(-age/tau)*(tau+age)))

print(norm_denom_real)
print(norm_denom, norm_denom_emma, norm_denom_new, norm_denom_wolf)


norm_emma = mass/tau*np.exp((-1)*age/tau)*(tau+age)+np.power(tau,2)

norm_lester = mass / ( (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2))

norm_real = mass / norm_denom_real

print(norm_emma, norm_lester, norm_real)

'''
('1A2744P', '1156_BEAGLE.fits.gz', 2255, 36, 0)
-35184372000000.0
(22932113000.0, 6485401.5)
'''

print('new')




age = 6485401.5
tau = 22932113000.0
mass = 1e9

time = np.linspace(0, 1e10, 1000)
sfr = time*np.exp(-time/tau)


plt.plot(time, sfr)
plt.plot((age, age), (0, max(sfr)))
plt.show()



sfr_function = lambda t: t*np.exp(-t/tau)
norm_denom_real = integrate.quad(sfr_function, 0, age)[0]




norm_denom_emma = tau*np.exp(-age/tau)*(tau+age)   +np.power(tau,2)

norm_denom = ( (-tau*np.exp(-age/tau)*(tau+age) )+np.power(tau,2))




print(norm_denom_emma, norm_denom, norm_denom_real)


norm_emma = mass/tau*np.exp((-1)*age/tau)*(tau+age)+np.power(tau,2)

norm_lester = mass / ( (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2))

norm_real = mass / norm_denom_real

print(norm_emma, norm_lester, norm_real)





# (-2199023300000.0, 5925514000.0, 1352515.1)

tau = 5925514000.0
age = 1352515.1

norm_denom = ( (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2))

print(norm_denom)

print(min(norm_denom))

for i in range(len(norm_denom)):
    if norm_denom[i] < 0:
        print(norm_denom[i], tau[i], age[i])


'''
HMMMM
(6189594600.0, 1670916.5)
(3.83110816094497e+19, 3.831108e+19)
(-6189594600.0, -6189594600.0)
(0.9997301, 0.9997301)
(6191266000.0, 6191266000.0)
(-3.8311084e+19, -3.8311084e+19)
(-2481614487552.0, -4398046500000.0)
'''


tau = 6189594600.0
age = 1670916.5

norm_denom = ( (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2))

print(norm_denom)



#%%

tau = 10**10.5 # 7, 10.5
age = 1e6 # 6, 10
mass = 1e5 # 5, 12

# 7 6 5 fine
# 10.5 6 5 fine
# 7 10 5


norm_denom2 = (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2)



sfr_function = lambda t: t*np.exp(-t/tau)

integral = integrate.quad(sfr_function, 0, age)
norm_denom = integral[0]
if integral[1]/integral[0] > 1e-5:
    print('integration error', i, file, integral[0], integral[1]/integral[0])
        




print(norm_denom, norm_denom2)




#%%

'''


log10_age = np.log10(age)
log10_mass = data['POSTERIOR PDF'].data['mass']            
log10_norm = log10_mass - np.log10(norm_denom)

exp_term = -age/tau
exp_idx = (exp_term <= -100.0)
exp_term[exp_idx] = -100.0
#            
temp_sfr = log10_norm + log10_age + np.log10(np.exp(exp_term))
sfr = np.where(temp_sfr < -30.0, -30.0, temp_sfr)

'''
























































