import sys
import numpy as np
import math

sys.path.append("/Users/user/code/pythonlib/mcurtislake-astrolibpy/astrolib/")
#import lumdist as ld

#works for flat UV slope!
def absMagCalc(appMag, z, lumD = None):
	if lumD is None:
		if hasattr(z, '__len__'):
			lumD = np.zeros(len(z))
			for i in range(len(z)):
				lumD[i] = ld.lumdist(z[i], silent=True)
		else:
			lumD = ld.lumdist(z, silent=True)
	return appMag - 5*np.log10(lumD*1.e6) + 5. + 2.5*np.log10(1+z)

def appMagToFnu(appMag):
	return 10**((appMag+48.6)/(-2.5))

def fnuToAppMag(fnu):
	return (-2.5)*np.log10(fnu)-48.6

def absMagToAppMag(absMag, z, lumD = None):
	if lumD is None:
		if hasattr(z, '__len__'):
			lumD = np.zeros(len(z))
			for i in range(len(z)):
				lumD[i] = ld.lumdist(z[i], silent=True)
		else:
			lumD = ld.lumdist(z, silent=True)
	return absMag + 5*np.log10(lumD*1.e6) - 5. - 2.5*np.log10(1+z)

def absMagToLnu(absMag, z):

#    if lumD is None:
#        if hasattr(z, '__len__'):
#            lumD = np.zeros(len(z))
#            for i in range(len(z)):
#                lumD[i] = ld.lumdist(z[i], silent=True)
#            else:
#                lumD = ld.lumdist(z, silent=True)
#    appMag = absMagToAppMag(absMag, z, lumD=lumD)
    fnu = appMagToFnu(absMag)

    #Lnu requires lumD in cm
    pcInCm = 3.08567758E18
    logLnu = np.log10(fnu * 4 * math.pi) +  2*np.log10(10. * pcInCm)

    return np.power(10,logLnu)

def appMagToLnu(appMag, z, lumD = None):

	#lumD calculated this way is in Mpc
	if lumD is None:
		if hasattr(z, '__len__'):
			lumD = np.zeros(len(z))
			for i in range(len(z)):
				lumD[i] = ld.lumdist(z[i], silent=True)
		else:
			lumD = ld.lumdist(z, silent=True)

	fnu = appMagToFnu(appMag)

	#Lnu requires lumD in cm
	pcInCm = 3.08567758E18
	Lnu = fnu * 4 * math.pi * (lumD * 1.E6 * pcInCm)**2

	return Lnu
 
 
	
