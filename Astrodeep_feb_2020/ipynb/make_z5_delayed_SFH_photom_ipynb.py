# =============================================================================
# 
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def integrate_delayed_history(scale, tau, age1, age2):
    mass = (-1)*scale*tau*np.exp((-1)*age2/tau)*(tau+age2) - \
           (-1)*scale*tau*np.exp((-1)*age1/tau)*(tau+age1)
    return mass

def find_age(mass, sfr, tau, minAge, maxAge):
    minScale = sfr/(minAge*np.exp((-1)*minAge/tau))
    minMass = integrate_delayed_history(minScale, tau, 0, minAge)
    maxScale = sfr/(maxAge*np.exp((-1)*maxAge/tau))
    maxMass = integrate_delayed_history(maxScale, tau, 0, maxAge)
#     print minMass, maxMass, mass
    if mass < minMass or mass > maxMass:
        #need to find new value of tau
        print 'here'
        return -1
    nAge = 500
    dAge = (maxAge-minAge)/nAge
    ageArr = np.fromiter((x for x in np.arange(minAge,maxAge+dAge,dAge)),np.float)

    scaleArr = sfr/(ageArr*np.exp((-1)*ageArr/tau))
    massArr = integrate_delayed_history(scaleArr, tau, 0, ageArr)

#     print ageArr[0], massArr[0], scaleArr[0], ageArr[-1], massArr[-1], scaleArr[-1]
#     print minAge, minMass, minScale, maxAge, maxMass, maxScale

# =============================================================================
#     plt.scatter(massArr, ageArr)
#     plt.show()
# =============================================================================
    #interpolate age vs. mass array to find age required to produce required mass
    #pylab.figure()
    #pylab.plot(np.log10(massArr), log10(ageArr))
    f = interpolate.interp1d(massArr, ageArr)
#     print '**', np.max(massArr), np.min(massArr), mass
#     print '***', np.max(ageArr), minAge, np.min(ageArr), maxAge
    age = f(mass)
    #pylab.title(str(log10(tau))+str(log10(mass))+" "+str(log10(age)))
    return age  


# =============================================================================
# 
# =============================================================================

nObj = 1000

minMass = 7.
maxMass = 12.

#assign masses drawn from broad Guassian distribution
Marr = np.random.normal(size=nObj,scale=0.5,loc=(minMass+maxMass)/2)

#rechoose a mass if not within minMass and maxMass
while np.max(Marr) > maxMass or np.min(Marr) < minMass:
    tempIdx = np.where((Marr > maxMass) | (Marr < minMass))[0]
    print(tempIdx)
    Marr[tempIdx] = np.random.normal(size=len(tempIdx),scale=0.5,loc=(minMass+maxMass)/2)
plt.figure()
plt.hist(Marr, bins=50)
plt.show()

#calculate sfr from eqn of straight line + scatter
intercept = -8.
slope = 1.
scatter = 0.3

Sfr = slope*Marr+intercept+np.random.normal(size=nObj,scale=scatter)
plt.scatter(Marr, Sfr)













import cosmolopy.distance as cd
import cosmolopy.constants as cc
cosmo = {'omega_M_0' : 0.286, 'omega_lambda_0' : 0.714, 'h' : 0.696}
cosmo = cd.set_omega_k_0(cosmo)
zObs=2.0
ageUniv = cd.age(zObs, **cosmo)/cc.yr_s

minAge = 1.E6
maxAge = ageUniv



from astropy.table import Table




minTau = 7.
maxTau = 10.5






tauArr = np.zeros(nObj)
scaleArr = np.zeros(nObj)
ageArr = np.zeros(nObj)

for j in range(nObj):
    #assign tau
    s = 10**Sfr[j]
    m = 10**Marr[j]
    tau = 10**(np.random.uniform(size=1)*(maxTau-minTau)+minTau)
    age = find_age(m, s, tau, minAge, maxAge)
    while age == -1:
        tau = 10**(np.random.uniform(size=1)*(maxTau-minTau)+minTau)
        age = find_age(m, s, tau, minAge, maxAge)
    tauArr[j] = np.log10(tau)
    ageArr[j] = np.log10(age)
    t = 10**ageArr[j]
    scaleArr[j] = s/(t*np.exp((-1)*t/tau))
    
    
    
outputDict = {}
outputDict['mass']=Marr
#print tauArr
outputDict['tau']=tauArr
outputDict['max_stellar_age']=ageArr
outputTable = Table(outputDict)
#outputTable.write("delayed_BEAGLEinput.fits", overwrite=True)


















