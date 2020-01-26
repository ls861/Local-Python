
import numpy as np
import matplotlib.pyplot as plt
import cosmolopy.distance as cd
import cosmolopy.constants as cc
import matplotlib.colors as mcolors

cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)

size        = np.int(1E5)

tau         = np.random.uniform(low=7,    high=10,    size=size)    # log
z           = np.random.uniform(low=4.5,    high=5.5,   size=size)  # redshift

age_h       = np.log10((10**9)*cd.age(4.5, **cosmo)/cc.Gyr_s)       # log, 9.107
age_z_15    = cd.age(15, **cosmo)/cc.Gyr_s                          # Gyr, age of universe at z=15
age         = np.random.uniform(low=6,      high=age_h, size=size)  # log, age of galaxy
age_Gyr     = (10**age)*(10**-9)                                    # Gyr, age of galaxy

t           = cd.age(z, **cosmo)/cc.Gyr_s                           # Gyr, same as Emma's age of universe
t0          = t-age_Gyr                                             # Gyr

mass    = np.random.uniform(low=5,      high=12,    size=size)      # log, mass of galaxy
massArr = np.arange(5, 12, 1)

#work out the sfr for all the sampled points

### TopHat ###

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
        elif t_yr - t0_yr > tau_yr:
            A = m / (t_yr - t0_yr)
        else:
            A = m / tau_yr

        sfr1 = A * np.heaviside(t_yr-t0_yr, 1) * (1 - np.heaviside((t_yr-t0_yr-tau_yr), 1))

        if sfr1 > 0:

            mass_used.append(mass[i])
            sfr.append(np.log10(sfr1))
        
plt.figure(figsize=(14, 10))
plt.hist2d(mass_used, sfr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")
plt.show()



'''
### Exponential ###

sfr = []
mass_used = []

for i in range(size):

    if t[i] - t0[i] < t[i] - age_z_15:
    
        m = 10**mass[i]             # solar masses
        t_yr = t[i]*(10**9)         # yrs
        t0_yr = t0[i]*(10**9)       # yrs
        tau_yr = 10**tau[i]         # yrs
        


        mass_used.append(mass[i])
        sfr.append(np.log10(A * np.heaviside(t_yr-t0_yr, 1) * np.exp(-(t_yr-t0_yr)/tau_yr)))
        
plt.figure(figsize=(14, 10))
plt.hist2d(mass_used, sfr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")
plt.show()


'''


'''
plt.figure(figsize=(14, 10))
plt.hist2d(mass_used, sfr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")


ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s
ageUniv_15 = 0.269

plt.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))
plt.plot(massArr, np.log10(10**massArr/(1E6)))


#plt.xlim(6,11)
#plt.ylim(-4,3)

plt.xlabel("log mass")
plt.ylabel("log sfr")
'''
#pylab.savefig("delayed_prior_z5_lester.pdf")
