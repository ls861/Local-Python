
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

age_z_4p5   = cd.age(4.5, **cosmo)/cc.yr_s                                      # yr, (log, 9.107)
age_z_15    = cd.age(15, **cosmo)/cc.yr_s                                       # yr, age of universe at z=15
age_galaxy  = (10**np.random.uniform(low=6, high=10, size=size))                # yr, age of galaxy

t_arr       = cd.age(z, **cosmo)/cc.yr_s                                        # yr, age of Universe
t0_arr      = t_arr-age_galaxy                                                  # yr, start of star formation
tau_arr     = np.random.uniform(low=6, high=11, size=size)                      # log, width of function

m_arr       = np.random.uniform(low=5, high=12, size=size)                      # log, mass of galaxy
massArr     = np.arange(5, 13, 1)

### Lognormal ###

sfr_arr     = []
mass_used   = []
int_err     = []

for i in range(size):

    if age_z_15 < t0_arr[i]:                    # age of galaxy is less than time available for star formation

        m       = 10**m_arr[i]                                                  # solar masses
        t       = t_arr[i]                                                      # yrs
        t0      = t0_arr[i]                                                     # yrs
        tau     = 10**tau_arr[i]                                                # yrs
        
        integrand = lambda T: np.heaviside(T-t0, 1) * (1/T) * np.exp( -(((np.log(T-t0))**2) / (2*(tau**2))) )
        integral  = quad(integrand, t0, t)
        
        A = m / integral[0]                     # solar masses / yr
        int_err.append(integral[1]/integral[0])
        
        sfr = A * np.heaviside(t-t0, 1) * (1/t) * np.exp( -(((np.log(t-t0))**2) / (2*(tau**2))) )

        mass_used.append(m_arr[i])                                              # log, mass of galaxy
        sfr_arr.append(np.log10(sfr))                                           # log, SFR

print(len(t_arr), len(sfr_arr), max(int_err))
        
plt.figure(figsize=(7, 5))
plt.hist2d(mass_used, sfr_arr, bins=[50,100], cmap='Blues', normed=True)
c=plt.colorbar()
c.set_label("weighting in prior")
plt.xlabel("log mass")
plt.ylabel("log sfr")

#plt.xlim(6,11)
#plt.ylim(-4,3)

ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.yr_s                                      # yr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.yr_s                                      # yr

plt.plot(massArr, np.log10(10**massArr/(1E6)), label='1E6')
plt.plot(massArr, np.log10(10**massArr/ageUniv_4p5), label='z=4.5')
plt.plot(massArr, np.log10(10**massArr/ageUniv_5p5), label='z=5.5')


plt.legend()
plt.show()

#pylab.savefig("delayed_prior_z5_lester.pdf")

