import sys
sys.path.append('/Users/efcl/pythonLib/usefulRoutines/')
import absMagCalc as amc
import numpy as np
import math
import copy
sys.path.append("/Users/efcl/pythonlib/mcurtislake-astrolibpy/astrolib/")
#import lumdist as ld

def conversion_factor_const(age):
    a=[-27.822, -0.969, 2.511, -2.0, 0.704, -0.089]
    cv = 0.
    for b in range(len(a)):
        cv = cv + a[b]*(1./(np.log(age/1.E6)))**b
    print 'age 2: ', age

    return 10**cv

#def conversion_factor_const(age):
#    a=np.array([2.32452038E-3, -5.80329863E-2, 3.18126368E-1, 2.7915886, -3.459012E1, 6.79967609E1])
#    f = np.poly1d(a)
#    cv = f(np.log10(age))
##    for b in range(len(a)):
##        cv = cv + a[b]*(np.log(age))**b
#
#    return 10**cv


def Reddy_UV_to_SFR(appMag, z, age, auv = 0, chab = True):
    print 'Reddy_UV_to_SFR AUV input: ', auv
    Lnu = amc.appMagToLnu(appMag, z)
    cv = conversion_factor_const(age)
    if chab == False:
      imf_conv = 1.
    else:
      imf_conv = 0.63
    sfr = cv * Lnu * 10**(0.4*auv)*imf_conv#/(1+z)
    return sfr

def Reddy_UV_to_SFR_abs(absMag, z, age, auv = 0, chab = True, lumD=None):
    print '** Reddy_UV_to_SFR AUV input: ', auv
    Lnu = amc.absMagToLnu(absMag, z)
    print 'age: ', age
    cv = conversion_factor_const(age)
    if chab == False:
        imf_conv = 1.
    else:
        imf_conv = 0.63
    sfr = cv * Lnu * 10**(0.4*auv)*imf_conv
    return sfr

def Kennicutt_UV_to_SFR(appMag, z, chab = True, absMag = False):
    if absMag:
        Lnu = amc.absMagToLnu(appMag, z)
    else:
        Lnu = amc.appMagToLnu(appMag, z)
    sfr = 1.4E-28 * Lnu/(1+z)
    if chab:
        sfr = sfr*0.63
    return sfr

def CB17_UV_to_SFR(appMag, z, chab = True, absMag = False, June16=True):
    if absMag:
        Lnu = amc.absMagToLnu(appMag, z)
    else:
        Lnu = amc.appMagToLnu(appMag, z)

    if June16:
        correctionFactor = 0.923571
    else:
        correctionFactor = 0.908679
    sfr = 1.4E-28*0.908679 * Lnu/(1+z)
#    if chab:
#        sfr = sfr*0.63
    return sfr


def Kennicutt_Ha_to_SFR(flux, z, chab = True, lumD = None, solar_luminosity = False):
    
    #lumD calculated this way is in Mpc
    if lumD is None:
        if hasattr(z, '__len__'):
            lumD = np.zeros(len(z))
            for i in range(len(z)):
                lumD[i] = ld.lumdist(z[i], silent=True)
            else:
                lumD = ld.lumdist(z, silent=True)

    if not solar_luminosity:
        pcInCm = 3.08567758E18
        logLum = np.log10(flux * 4 * math.pi) + 2*np.log10(lumD * 1.E6 * pcInCm)
    else:
        logLum = np.log10(flux) + np.log10(3.826e33)
    sfr = np.power(10,np.log10(7.9E-42)+logLum)
    if chab:
        print 'introducing chab conversion in ha SFR'
        sfr = sfr*0.63
    return sfr

def CB17_Ha_to_SFR(flux, z, chab = True, lumD = None, solar_luminosity = False, June16 = True):
    #NOTE TO SELF: Not sure this is correct, should re-do and verify conversion factors as I think
    #They already include the difference between salpeter and chabrier!
    
    #lumD calculated this way is in Mpc
    if lumD is None:
        if hasattr(z, '__len__'):
            lumD = np.zeros(len(z))
            for i in range(len(z)):
                lumD[i] = ld.lumdist(z[i], silent=True)
            else:
                lumD = ld.lumdist(z, silent=True)

    if not solar_luminosity:
        pcInCm = 3.08567758E18
        logLum = np.log10(flux * 4 * math.pi) + 2*np.log10(lumD * 1.E6 * pcInCm)
    else:
        logLum = np.log10(flux) + np.log10(3.826e33)
    if June16:
        correctionFactor = 0.6906
    else:
        correctionFactor = 0.696264
    sfr = np.power(10,np.log10(7.9E-42*correctionFactor)+logLum)
#    if chab:
#        print 'introducing chab conversion in ha SFR'
#        sfr = sfr*0.63
    return sfr


def MeurerCorrection(beta, uvMag, redshift, flambda=False):
    if flambda:
        uvFlambda = np.array(uvMag)*(1600/1500)**beta
        uvFnu = uvFlambda*((1600.*(1+redshift))**2/2.998E18)
        uvMag = (-2.5)*np.log10(uvFnu) - 48.58
    A1600 = 4.43 + 1.99*beta
    tempIdx = np.where(beta < 0)
    uvCorr = uvMag
    uvCorr[tempIdx] = uvMag[tempIdx] - A1600[tempIdx]
    if not flambda:
        return uvCorr
    else:
        print uvFlambda[0:10], uvCorr[0:10]
        uvCorrFnu = 10**((uvMag + 48.58)/(-2.5))
        print uvCorrFnu[0:10]
        uvCorrFlambda = uvCorrFnu/((1500.*(1+redshift))**2/2.998E18)
        test = uvFnu/10**(A1600/(-2.5))
        test2 = test/((1500.*(1+redshift))**2/2.998E18)
        print beta[0:10], test2[0:10], uvCorrFlambda[0:10]
        return uvCorrFlambda

def caseBcorrection(Ha, Hb):
    EBV = (2.5*np.log10(Ha/Hb)-2.5*np.log10(2.86))/1.07
    Aha = 2.45*EBV
    temp = Ha/10**(Aha/(-2.5))
    idxArr = np.all([Ha > 0, Hb > 0], axis=0)
    HaCorr = copy.deepcopy(Ha)
    HaCorr[idxArr] = temp[idxArr]
    return HaCorr

def Santini_SFR(appMagArr, appMagErr, lambdaArr, z):
    c=2.998E18 #A/s
    appMagArr = np.array(appMagArr)
    appMagErr = np.array(appMagErr)
    lambdaArr = np.array(lambdaArr)
    #first fit a straight line through the points M = -2.5(beta+2)*log(lambda)+c
    results = np.polyfit(np.log10(lambdaArr), appMagArr, 1)
    beta = results[0]/(-2.5)-2
    #first calculate the corrected luminosity at 1600A
    obs_1600 = 1600*(1+z)
    mag_1600 = results[1]+results[0]*np.log10(obs_1600)
    #correct luminosity using Meurer 1999 A1600 = 4.43+1.99*beta
    Meurer_correction = 4.43 + 1.99*beta
    Lnu = amc.appMagToLnu(mag_1600, z)
    Lnu_corr = amc.appMagToLnu(mag_1600-Meurer_correction, z)
    #Then use assumed intrinsic UV slope of -2.23 assumed for the Meurer relation to find Lnu_corr_1500
    beta_int = -2.23
    Lnu_corr_1500 = Lnu_corr*(1500**(beta_int-2))/(1600**(beta_int-2))
    #Kennicutt 2012 conversion factor log(SFR) = log(L)-log(cf)
    log_cf = 43.35
    log_sfr = np.log10(c/1500*Lnu_corr_1500/(1+z))-log_cf
    temp_sfr = np.log10(c/1500*Lnu/(1+z))-log_cf
    #use error in flux closest to rest-frame 1600 for error on M_1600 (bit of a cheat)
    delta = np.abs(lambdaArr-obs_1600)
    minIdx = np.argmin(delta)
    magErr = appMagErr[minIdx]
    Lnu_err = amc.appMagErrToLnuErr(appMagArr[minIdx],appMagErr[minIdx], z)
    log_Lnu_err = (0.434*(Lnu_err/Lnu))
    log_sfr_err = np.sqrt(log_Lnu_err**2+0.55**2)#adding in scatter in Meurer relation in quadrature
    return log_sfr, log_sfr_err, Meurer_correction

def Santini_SFR_no_err(appMagArr, appMagErr, lambdaArr, z):
    c=2.998E18 #A/s
    appMagArr = np.array(appMagArr)
    appMagErr = np.array(appMagErr)
    lambdaArr = np.array(lambdaArr)
    #first fit a straight line through the points M = -2.5(beta+2)*log(lambda)+c
    results = np.polyfit(np.log10(lambdaArr), appMagArr, 1)
    beta = results[0]/(-2.5)-2
    #first calculate the corrected luminosity at 1600A
    obs_1600 = 1600*(1+z)
    mag_1600 = results[1]+results[0]*np.log10(obs_1600)
    #correct luminosity using Meurer 1999 A1600 = 4.43+1.99*beta
    Meurer_correction = 4.43 + 1.99*beta
    Lnu = amc.appMagToLnu(mag_1600, z)
    Lnu_corr = amc.appMagToLnu(mag_1600-Meurer_correction, z)
    #Then use assumed intrinsic UV slope of -2.23 assumed for the Meurer relation to find Lnu_corr_1500
    beta_int = -2.23
    Lnu_corr_1500 = Lnu_corr*(1500**(beta_int-2))/(1600**(beta_int-2))
    #Kennicutt 2012 conversion factor log(SFR) = log(L)-log(cf)
    log_cf = 43.35
    log_sfr = np.log10(c/1500*Lnu_corr_1500/(1+z))-log_cf
    temp_sfr = np.log10(c/1500*Lnu/(1+z))-log_cf
    #use error in flux closest to rest-frame 1600 for error on M_1600 (bit of a cheat)
    delta = np.abs(lambdaArr-obs_1600)
    minIdx = np.argmin(delta)
    magErr = appMagErr[minIdx]
    return log_sfr
