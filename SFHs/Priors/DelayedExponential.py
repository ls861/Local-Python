
import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
#import matplotlib.colors as mcolors

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)
z           = np.random.uniform(low=4.5,    high=5.5,   size=size)          # redshift

age_z_4p5_Gyr = cd.age(4.5, **cosmo)/cc.Gyr_s                               # Gyr
age_z_4p5   = np.log10((10**9)*age_z_4p5_Gyr)                               # log, 9.107

age_z_15    = cd.age(15, **cosmo)/cc.Gyr_s                                  # Gyr, age of universe at z=15

age         = np.random.uniform(low=6, high=np.log10((10**9)*(age_z_4p5_Gyr-age_z_15)), size=size)      # log, age of galaxy
age_Gyr     = (10**age)*(10**-9)                                            # Gyr, age of galaxy

t           = cd.age(z, **cosmo)/cc.Gyr_s                                   # Gyr, same as Emma's age of universe
t0          = t-age_Gyr                                                     # Gyr, start of star formation
tau         = np.random.uniform(low=8,      high=12,    size=size)          # log, width of tophat

mass        = np.random.uniform(low=5,      high=12,    size=size)          # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

#work out the sfr for all the sampled points

### Delayed Exponential ###


print(age_z_4p5)
print(np.log10((10**9)*age_z_15))

sfr = []
mass_used = []

for i in range(size):

    if t[i] - t0[i] < t[i] - age_z_15:
    
        m = 10**mass[i]             # solar masses
        t_yr = t[i]*(10**9)         # yrs
        t0_yr = t0[i]*(10**9)       # yrs
        tau_yr = 10**tau[i]         # yrs
        
        if t_yr - t0_yr < 0:
            A = 0
            
        else:
            A = m / (-tau_yr* ( ((np.exp(-(t_yr-t0_yr)/tau_yr))*(t_yr-t0_yr+tau_yr)) -tau_yr) )

        sfr1 = A * np.heaviside(t_yr-t0_yr, 1) * (t_yr-t0_yr) * np.exp(-(t_yr-t0_yr)/tau_yr)

        if sfr1 > 0:

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

