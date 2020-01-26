
import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.special import erf
#import matplotlib.colors as mcolors

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)
z           = np.random.uniform(low=4.5,    high=5.5,   size=size)          # redshift

age_z_4p5   = np.log10((10**9)*cd.age(4.5, **cosmo)/cc.Gyr_s)               # log, 9.107
age_z_15    = cd.age(15, **cosmo)/cc.Gyr_s                                  # Gyr, age of universe at z=15

age_Gyr     = (10**np.random.uniform(low=6, high=age_z_4p5, size=size))*(10**-9)    # Gyr, age of galaxy

t           = cd.age(z, **cosmo)/cc.Gyr_s                                   # Gyr, same as Emma's age of universe
t0          = t-age_Gyr                                                     # Gyr, start of star formation
tau         = np.random.uniform(low=10,     high=11,    size=size)          # log, width of guassian
tpeak       = np.random.uniform(low=4,      high=10,    size=size)          # log, peak of guassian

mass        = np.random.uniform(low=5,      high=12,   size=size)           # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

#work out the sfr for all the sampled points

### Guassian ###

sfr = []
mass_used = []


test = 0
for i in range(size):

    if t[i] - t0[i] < t[i] - age_z_15:

        m = 10**mass[i]             # solar masses
        t_yr = t[i]*(10**9)         # yrs
        t0_yr = t0[i]*(10**9)       # yrs
        tpeak_yr = 10**tpeak[i]     # yrs
        tau_yr = 10**tau[i]         # yrs
        
        A = m / (-tau_yr*((np.pi/2)**0.5)*((erf((tpeak_yr-t_yr)/((2**0.5)*tau_yr)))-(erf((tpeak_yr-t0_yr)/((2**0.5)*tau_yr)))))
        
        sfr1 = A * np.exp( -0.5*(((t_yr-tpeak_yr) / tau_yr)**2) )

        if sfr1 > 0:
            test+=1
            #print(test)
            mass_used.append(mass[i])
            sfr.append(np.log10(sfr1))
          
print(len(t))
print(len(sfr))
        
plt.figure(figsize=(14, 10))
plt.hist2d(mass_used, sfr, bins=[50,100], cmap='Blues', normed=True)
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

