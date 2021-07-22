import numpy as np
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt

'''
args='-r /Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/ --santini --input-cat /Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/data/astrodeep_A2744_p_subset_RF1_001.fits'
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
    return log_sfr, log_sfr_err, Meurer_correction

def appMagToLnu(appMag, z):
    import math
	#lumD calculated this way is in Mpc
    lumD = lumdist(z)
    fnu = appMagToFnu(appMag)

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
    return lumD.value

def run_main(args):
  if args.santini:
    inputCat = fits.open(args.inputCat)
    summCat = args.resultsDir+"pyp-beagle/data/BEAGLE_summary_catalogue.fits"
    data_fits = fits.open(summCat)
    idArr = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    


    ### UPDATE TO MATCH SANTINI EVEN MORE ###
#    redshift_med = np.asarray(data_fits['GALAXY PROPERTIES'].data['redshift_median'], dtype=np.float64)
    redshift_med = inputCat[1].data['ZBEST'][idArr-1]

    if args.Mtot:
        mass_med = np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_tot_median'], dtype=np.float64))
    else:
        mass_med = np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_star_median'], dtype=np.float64))

#    temp_sfr = np.asarray(data_fits['STAR FORMATION'].data['SFR_median'], dtype=np.float64)
#    data['sfr'] = {'value':np.log10( np.where(temp_sfr==0, 1e-30, temp_sfr) )}
    data_fits.close()
    
    SFR_santini_ID = []
    SFR_santini = []
    SFR_err_santini = []
    Meurer_correction_arr = []
    #find filters spanning 1280-2600A
    
    filter_label = np.array(['b_B435', 'b_V606', 'b_I814', 'b_Y105', 'b_J125', 'b_JH140', 'b_H160', 'b_Ks', 'b_CH1', 'b_CH2'])
    filter_fwhm_centre = np.array([4348.65, 5926.47, 7975.65, 10530.87, 12495.71, 13976.13, 15433.07, 21440.35, 35465.62, 45024.31])
    filter_fwhm = np.array([939, 2322.94, 1856, 2917.03, 3005.2, 3940.88, 2874.18, 3249.92, 7431.71, 10096.82])   
    filter_low = filter_fwhm_centre - filter_fwhm/2
    filter_high = filter_fwhm_centre + filter_fwhm/2
        
    
#    plt.hist(redshift_med, bins=30)
#    plt.show()
    print(len(redshift_med), len(redshift_med[abs(redshift_med-1.65)<0.35]))
    for i in range(len(mass_med)):
        tempIdx = (np.where(inputCat[1].data['ID'] == idArr[i]))[0][0]
#        print(tempIdx)
            
        
        
        
        lambda_1280 = 1280*(redshift_med[i]+1.0)
        lambda_2600 = 3100*(redshift_med[i]+1.0) # 3100 ensures 2 filters available at z==1.3
        
        
#        l1280 = 1280*(1.3+1.0)
#        l2600 = 3100*(1.3+1.0)
#        
#        h1280 = 1280*(2.0+1.0)
#        h2600 = 3100*(2.0+1.0)       
        

        temp_fnu_lester = []
        appMagArr = []
        appMagErr = []
        lambdaArr = []
        
#        print(filter_high)
#        print(filter_low)       
#        print(l1280, l2600, h1280, h2600)
        
        
#        print(filterInfo[1].header['TFIELDS'])
        

        for f in range(len(filter_label)):
            
            if lambda_1280  < filter_low[f] and lambda_2600 > filter_high[f]:
                
                temp_fnu = inputCat[1].data[filter_label[f]][tempIdx]
                
                if temp_fnu > 0:
                    
                    temp_fnu_err = inputCat[1].data[filter_label[f].replace('_', '_err')][tempIdx]
                    appMagArr.append(-2.5 * np.log10(temp_fnu*1e-6) + 8.90)
                    appMagErr.append(1.0857*temp_fnu_err/temp_fnu)
                    lambdaArr.append(filter_fwhm_centre[f])
        
                    temp_fnu_lester.append(temp_fnu)
#        if len(appMagArr) < 2:
#            print ("error: ", idArr[i], lambdaArr, redshift_med[i], appMagArr) 

#        if len(appMagArr) < 2:
#            sys.exit()

        if len(appMagArr) > 1:
            temp_sfr, temp_sfr_err, meurer_correction = Santini_SFR(appMagArr, appMagErr, lambdaArr, redshift_med[i])
            
#            temp_sfr, meurer_correction = Santini_SFR(appMagArr, appMagErr, lambdaArr, redshift_med[i])            
            SFR_santini_ID.append(idArr[i])
            SFR_santini.append(temp_sfr)
            SFR_err_santini.append(temp_sfr_err)
            Meurer_correction_arr.append(meurer_correction)

    SFR_santini = np.array(SFR_santini)
    SFR_err_santini = np.array(SFR_err_santini)


    idx_santini = np.isin(idArr, SFR_santini_ID)

    print(len(redshift_med[idx_santini]), len(redshift_med[idx_santini][abs(redshift_med[idx_santini]-1.65)<0.35]))
    
    
#    plt.scatter(mass_med[idx_santini], SFR_santini)
#    plt.show()
    

    print(np.polyfit(mass_med[idx_santini], SFR_santini, 1))
    
    print(len(SFR_santini))
    
    sbf = './results/{}_'.format(args.resultsDir[-2])
    np.save(sbf+'mass_santini.npy', mass_med[idx_santini])
    np.save(sbf+'sfr_santini.npy', SFR_santini)
    np.save(sbf+'sfr_err_santini.npy', SFR_err_santini)
    np.save(sbf+'sfr_ID_santini.npy', SFR_santini_ID)
    np.save(sbf+'meurer_correction.npy', Meurer_correction_arr)  
    
    



parser = argparse.ArgumentParser()

parser.add_argument(
    '-r', '--results-dir',
    help="Directory containing the fits.gz output files from BEAGLE fitting",
    action="store",
    type=str,
    dest="resultsDir",
    required=True
)

parser.add_argument(
    '--m-tot',
    help="Measure results from M_tot rather than M_star",
    action="store_true",
    dest="Mtot",
    required=False,
    default=False
)

parser.add_argument(
  '--santini',
  help="Measure results with Santini+17 method",
  action="store_true",
  dest="santini",
  required=False,
  default=False
)

parser.add_argument(
  '--input-cat',
  help="Input photometric catalogue, required if --santini flag used",
  action="store",
  type=str,
  dest="inputCat",
  required=True,
  default=None
)

args = parser.parse_args()

f = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/'
fields = ['A2744_c','A2744_p','M0416_c','M0416_p','M0717_c','M0717_p','M1149_c','M1149_p']
#fields = ['A2744_c']

for i in range(len(fields)):
    args.resultsDir = '{}fields/{}/'.format(f, str(i))
    args.inputCat = '{}data/astrodeep_{}_subset_RF1_001.fits'.format(f, fields[i])
    
#    print(args.resultsDir)
#    print(args.inputCat)
    
    run_main(args)
    
    



