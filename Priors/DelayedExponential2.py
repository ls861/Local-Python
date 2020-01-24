
import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.special import erf
from scipy.integrate import quad
import matplotlib.colors as mcolors


cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)
z           = np.random.uniform(low=4.5, high=5.5, size=size)                   # redshift

age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = 10**np.random.uniform(low=6, high=10, size=size)                  # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation
tau_arr     = 10**np.random.uniform(low=7.5, high=10, size=size)                # yr, width of function

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)


### Delayed Exponential ###

sfr_arr = []
mass_used = []

ssfr_arr = []
tau_log10_arr = []

for i in range(size):

    if age_z_15 < t0_arr[i]:
    
        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0      = t0_arr[i]                                                     # yrs
        tau     = tau_arr[i]                                                    # yrs

        A = m / (-tau* ( ((np.exp(-(t-t0)/tau))*(-t0+tau+t)) - tau) )
        # https://www.wolframalpha.com/input/?i=%28t-B%29*exp%28-%28t-B%29%2FC%29  
        
        sfr = A * (t-t0) * np.exp(-(t-t0)/tau)
        mass_used.append(m_arr[i])
        sfr_arr.append(np.log10(sfr))
        
        ssfr_arr.append(np.log10(sfr/m))
        tau_log10_arr.append(np.log10(tau))

        '''
        xax = (10**9)*np.arange(0.1, 15, 0.01)
        plt.figure(figsize=(10, 5))
        plt.plot(xax, A * (xax-t0) * np.exp(-(xax-t0)/tau))
        plt.plot((t, t), (0, 1E100))
        plt.plot((t0, t0), (0, 1E10))
        plt.xlim((0, 14E9))
        plt.ylim((0, max(A * (xax-t0) * np.exp(-(xax-t0)/tau))))
        plt.show()
        print(A, m)
        '''
            
print(len(t_arr))
print(len(sfr_arr))

### ### ### ### ###

plt.figure(figsize=(7, 5))
plt.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")

plt.xlim(6,11)
plt.ylim(-4,3)

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s

plt.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(1E6)))

plt.show()

### ### ### ### ###

plt.figure(figsize=(7, 5))
plt.hist2d(ssfr_arr, tau_log10_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log ssfr")
plt.ylabel("log tau")

#plt.xlim(6,11)
#plt.ylim(-4,3)

plt.show()

### ### ### ### ###

10**np.random.uniform(low=7.5, high=10, size=size)
tau_calc_arr = 10**np.arange(7.5, 10, 0.01)
mass_calc_arr = np.zeros(size)

A = 1
time_now = 1.12 * 1E9
t = time_now
msa = 5E8
t0 = t - msa

mass_calc_arr = A*(-tau_calc_arr* ( ((np.exp(-(t-t0)/tau_calc_arr))*(-t0+tau_calc_arr+t)) - tau_calc_arr) )


plt.figure(figsize=(10, 5))
plt.plot(tau_calc_arr, mass_calc_arr)
time_now = 1.12 * 1E9
plt.plot((msa, msa), (0, 13E16), color='k', linestyle=':')   
plt.xlim(0, 1E10)
plt.ylim(0, 13E16)
plt.xlabel('TAU')
plt.ylabel('MASS')
plt.show()

plt.figure(figsize=(10, 5))

time = 1E9*np.linspace(0, 10, 1E6)

plt.plot(time, A * (time-t0) * np.exp(-(time-t0)/msa))
plt.plot((time_now, time_now), (0, 2E8), color='k', linestyle=':')   
plt.xlim(0, 1E10)
plt.ylim(0, 2E8)
plt.xlabel('AGE OF UNIVERSE')
plt.ylabel('SFR')
plt.show()





### Delayed Exponential ###
'''

print(np.log10((10**9)*age_z_15))

sfr_arr = []
mass_used = []

for i in range(size):

    if age_z_15 < t0_arr[i]:
    
        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0      = t0_arr[i]                                                     # yrs
        tau     = tau_arr[i]                                                    # yrs

        if t - t0 < 0:
            A = 0
        else:
            A = m / (-tau* ( ((np.exp(-(t-t0)/tau))*(t-t0+tau)) - tau) )

        sfr = A * np.heaviside(t-t0, 1) * (t-t0) * np.exp(-(t-t0)/tau)

        if sfr > 0:

            mass_used.append(m_arr[i])
            sfr_arr.append(np.log10(sfr))

'''

'''    
integrand = lambda T: 1 / (((T/tau)**alpha)+((T/tau)**-beta))
integral  = quad(integrand, 0, t)

A = m / integral[0]                                                 # solar masses / yr
#int_err.append(integral[0]/integral[1])
sfr = A / (((t/tau)**alpha)+((t/tau)**-beta))
mass_used.append(m_arr[i])                                          # log, mass of galaxy
sfr_arr.append(np.log10(sfr))             


plt.plot(np.arange(0.1, 15, 0.01), A / (((((10**9)*np.arange(0.1, 15, 0.01))/tau)**alpha)+((((10**9)*np.arange(0.1, 15, 0.01))/tau)**-beta)))
plt.plot((t*1E-9, t*1E-9), (0, A))
plt.plot((t0*1E-9, t0*1E-9), (0, A))
plt.xlim(0,2)
plt.show()
print(m)
'''




'''
            
print(len(t_arr))
print(len(sfr_arr))
        
plt.figure(figsize=(7, 5))
plt.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")

#plt.xlim(6,11)
#plt.ylim(-4,3)

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s

plt.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(1E6)))

plt.show()

#pylab.savefig("delayed_prior_z5_lester.pdf")

'''