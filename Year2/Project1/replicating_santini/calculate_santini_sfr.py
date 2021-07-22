import numpy as np
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt

'''
UPDATE as this should only require the AD catalogue
should use ZBEST as redshift
'''

def Santini_SFR(appMagArr, appMagErr, lambdaArr, z):
    c=2.998E18 #A/s
    appMagArr = np.array(appMagArr)
    appMagErr = np.array(appMagErr)
    lambdaArr = np.array(lambdaArr)
    # first fit a straight line through the points M = -2.5(beta+2)*log(lambda)+c
    results = np.polyfit(np.log10(lambdaArr), appMagArr, 1)
    beta = results[0]/(-2.5)-2
    #first calculate the corrected luminosity at 1600A
    obs_1600 = 1600*(1+z)
    mag_1600 = results[1]+results[0]*np.log10(obs_1600)
    #correct luminosity using Meurer 1999 A1600 = 4.43+1.99*beta
    Meurer_correction = 4.43 + 1.99*beta
    Lnu = appMagToLnu(mag_1600, z)
    print(mag_1600, mag_1600-Meurer_correction)
    Lnu_corr = appMagToLnu(mag_1600-Meurer_correction, z)
    #Then use assumed intrinsic UV slope of -2.23 assumed for the Meurer relation to find Lnu_corr_1500
    beta_int = -2.23
    Lnu_corr_1500 = Lnu_corr*(1500**(beta_int-2))/(1600**(beta_int-2))
    #Kennicutt 2012 conversion factor log(SFR) = log(L)-log(cf)
    log_cf = 43.35
    log_sfr = np.log10(c/1500*Lnu_corr_1500/(1+z))-log_cf
#    temp_sfr = np.log10(c/1500*Lnu/(1+z))-log_cf
    #use error in flux closest to rest-frame 1600 for error on M_1600 (bit of a cheat)
    delta = np.abs(lambdaArr-obs_1600)
    minIdx = np.argmin(delta)
#    magErr = appMagErr[minIdx]
    Lnu_err = appMagErrToLnuErr(appMagArr[minIdx], appMagErr[minIdx], z)
    log_Lnu_err = (0.434*(Lnu_err/Lnu))
    log_sfr_err = np.sqrt(log_Lnu_err**2+0.55**2)#adding in scatter in Meurer relation in quadrature
    print(Lnu_corr)
    return log_sfr, log_sfr_err, Meurer_correction

def appMagToLnu(appMag, z):
    import math
	#lumD calculated this way is in Mpc
    lumD = lumdist(z)
    fnu = appMagToFnu(appMag)
#    print(fnu)
    #Lnu requires lumD in cm
    pcInCm = 3.08567758E18
    Lnu = fnu * 4 * math.pi * (lumD * 1.E6 * pcInCm)**2
    return Lnu

def appMagErrToLnuErr(appMag, appMagErr, z):
    import math
    #lumD calculated this way is in Mpc
    lumD = lumdist(z)
    fnu = appMagToFnu(appMag)
    fnuErr = appMagErr*fnu/1.0857
    
    #Lnu requires lumD in cm
    pcInCm = 3.08567758E18
    LnuErr = fnuErr * 4 * math.pi * (lumD * 1.E6 * pcInCm)**2
    return LnuErr
 
def appMagToFnu(appMag):
    return 10**((appMag+48.6)/(-2.5))

def lumdist(z):
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    lumD = cosmo.luminosity_distance(z)   
#    print(lumD)
    return lumD.value


AD_location = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/astrodeep_rawfile_1234_ABCZ.npy'
AD = np.load(AD_location)
#print(AD.dtype.names)
#print(len(AD))

# SUBSET FOR TESTING
AD = AD[:20]

redshift_med = AD['ZBEST']

SFR_santini = []
SFR_err_santini = []
Meurer_correction_arr = []

#find filters spanning 1280-2600A
filter_label = np.array(['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'])
filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])   
filter_low = filter_fwhm_centre - filter_fwhm/2
filter_high = filter_fwhm_centre + filter_fwhm/2

for i in range(len(redshift_med)):
    
    lambda_1280 = 1280*(redshift_med[i]+1.0)
    lambda_2600 = 3100*(redshift_med[i]+1.0) # 3100 ensures 2 filters available at z==1.3
    
    temp_fnu_lester = []
    appMagArr = []
    appMagErr = []
    lambdaArr = []
    
    for f in range(len(filter_label)):
        
        if lambda_1280  < filter_low[f] and lambda_2600 > filter_high[f]:
            
            temp_fnu = AD[filter_label[f]][i]
            
            if temp_fnu > 0:
                
                temp_fnu_err = AD[filter_label[f].replace('_', '_err')][i]
                appMagArr.append(-2.5 * np.log10(temp_fnu*1e-6) + 8.90)
                appMagErr.append(1.0857*temp_fnu_err/temp_fnu)
                lambdaArr.append(filter_fwhm_centre[f])

    if len(appMagArr) >= 2:
        temp_sfr, temp_sfr_err, meurer_correction = Santini_SFR(appMagArr, appMagErr, lambdaArr, redshift_med[i])
                  
        SFR_santini.append(temp_sfr)
        SFR_err_santini.append(temp_sfr_err)
        Meurer_correction_arr.append(meurer_correction)
#        print(temp_sfr)
    else:
        SFR_santini.append(-99.0)
        SFR_err_santini.append(-99.0)
        Meurer_correction_arr.append(-99.0)        

SFR_santini = np.array(SFR_santini)
SFR_err_santini = np.array(SFR_err_santini)
Meurer_correction_arr = np.array(Meurer_correction_arr)


'''
sbf = './calculate_santini_sfr_results/'
np.save(sbf+'sfr_santini.npy', SFR_santini)
np.save(sbf+'sfr_err_santini.npy', SFR_err_santini)
np.save(sbf+'meurer_correction.npy', Meurer_correction_arr)  
'''






