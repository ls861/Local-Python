import numpy as np
import pylab

#size is the number of points you're going to sample from the prior
#You effectively end up plotting a 2d histogram so making this number larger
#will make it look less noisy but at the expense of speed of running the
#script
size = np.int(1E5)

#set up the prior on tau (you could try plotting different priors, like
#the uniform prior on 1/tau (commented out below) rather than uniform prior on log tau which
#is set up here
minTau = 7.5 #these are log values
maxTau = 10
tau = np.random.uniform(size=size)*(maxTau-minTau)+minTau

#minTau = 1./10**10.5
#maxTau = 1./10**7.
#invTau = np.random.uniform(size=size)*(maxTau-minTau)+minTau
#tau = np.log10(1./invTau)

#let's set the redshift - you could make this a single value even
#The figures in my paper are plotted with z=5
minZ = 4.5
maxZ = 5.5
z = np.random.uniform(size=size)*(maxZ-minZ)+minZ
ageUniv = np.zeros(size)

import cosmolopy.distance as cd
import cosmolopy.constants as cc
cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
cosmo = cd.set_omega_k_0(cosmo)
for i in range(size):
    ageUniv[i] = cd.age(z[i], **cosmo)/cc.Gyr_s

#prior on log age
minAge = 6
maxAge = 12
age = np.random.uniform(size=size)*(maxAge-minAge)+minAge

#prior on log mass
minMass = 5
maxMass = 12
mass = np.random.uniform(size=size)*(maxMass-minMass)+minMass

#work out the sfr for all the sampled points
sfr = []
mass_used = []
age_z_15 = cd.age(15, **cosmo)/cc.Gyr_s
for i in range(size):
    if age[i] < np.log10((ageUniv[i] - age_z_15)*1E9):
        m = 10**mass[i]
        a = 10**age[i]
        t = 10**tau[i]
        mass_used.append(mass[i])
        scale = m/((-1)*t*np.exp((-1)*a/t)*(t+a)+t*np.exp(-0/t)*(t+0))
        sfr.append(np.log10(scale*a*np.exp((-1)*a/t)))
#print scale, sfr[-1], a*np.exp((-1)*a/t), t*np.exp((-1)*a/t)*(t+a)

import matplotlib.colors as mcolors
pylab.figure()
#pylab.hist2d(mass_used, sfr, bins=[50,50], norm=mcolors.LogNorm(), cmap='Blues')
pylab.hist2d(mass_used, sfr, bins=[50,100], cmap='Blues', normed=True)
c=pylab.colorbar()
c.set_label("weighting in prior")


ageUniv_5p5 = cd.age(5.5, **cosmo)/cc.Gyr_s #Gyr
ageUniv_4p5 = cd.age(4.5, **cosmo)/cc.Gyr_s
ageUniv_15 = 0.269

massArr = np.arange(minMass,maxMass,1)
pylab.plot(massArr, np.log10(10**massArr/(ageUniv_5p5*1E9)))
pylab.plot(massArr, np.log10(10**massArr/(ageUniv_4p5*1E9)))


pylab.plot(massArr, np.log10(10**massArr/(1E6)))
pylab.ylim(-4,3)
#pylab.xlim(6,11)
pylab.xlabel("log mass")
pylab.ylabel("log sfr")


pylab.savefig("delayed_prior_z5.pdf")
