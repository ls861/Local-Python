import numpy as np
import pylab
from astropy.io import fits
import os
import argparse
#from beagle_summary_catalogue import get1DInterval
import pickle
import numpy.polynomial.polynomial as poly
import sys
#sys.path.append('/nethome/curtis/mycode/pythonlib/usefulRoutines/')
#import Reddy_SFR_conversions_LS as SFR_conv
import copy
import matplotlib.pyplot as plt


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
    Lnu_err = appMagErrToLnuErr(appMagArr[minIdx],appMagErr[minIdx], z)
    log_Lnu_err = (0.434*(Lnu_err/Lnu))
    log_sfr_err = np.sqrt(log_Lnu_err**2+0.55**2)#adding in scatter in Meurer relation in quadrature
    return log_sfr, log_sfr_err, Meurer_correction
#    return log_sfr, Meurer_correction



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

  if args.salmon or args.santini:
    if args.inputCat == None:
      print 'inputCat must be specified if trying to use the Salmon+2014 measurement'
      return
    else:
      if args.filterInfo == None:
        print 'filterInfo must be specified if trying to use the Salmon+2014 measurement'
        return
      else:
        inputCat = fits.open(args.inputCat)
        filterInfo = fits.open(args.filterInfo)

  '''
  fileTest = os.path.isfile(args.resultsDir+"medianDict.p")
  if fileTest:
    outputDict = pickle.load(open(args.resultsDir+"medianDict.p","r"))
    idArr = outputDict['idArr']
    mass_med = outputDict['mass_med'] 
    mass_err = outputDict['mass_err']
    SFR_med = outputDict['SFR_med']
    SFR_err = outputDict['SFR_err']
    AUV_med = outputDict['AUV_med']
    AUV_err = outputDict['AUV_err']
    redshift_med = outputDict['redshift_med']
    redshift_err = outputDict['redshift_err']
    age_med = outputDict['age_med']
    age_err = outputDict['age_err']
  else: 
    fileList = os.listdir(args.resultsDir)

    idArr = []
    mass_med = []
    mass_err = []
    SFR_med = []
    SFR_err = []
    AUV_med = []
    AUV_err = []
    redshift_med = []
    redshift_err = []
    age_med = []
    age_err = []

    levels = [68]
    for file in fileList:
      if '.fits.gz' in file:
        print file
        idArr.append(np.int(file.replace("_BEAGLE.fits.gz","")))
        print idArr[-1]
        results = fits.open(args.resultsDir+file)
        prob = results['POSTERIOR PDF'].data['probability']
        
        
        if (args.Mstar):
          temp_arr = np.log10(results['GALAXY PROPERTIES'].data['M_star'])
        else:
          temp_arr = results['POSTERIOR PDF'].data['mass']
                
        mean, median, interval = get1DInterval(temp_arr, prob, levels)
        mass_med.append(median)
        mass_err.append(interval[0])

        if args.delayed:
          #instantaneous SFR if delayed
          age = np.float64(data['STAR FORMATION'].data['max_stellar_age'])
      	  tau = np.float64(data['STAR FORMATION BINS'].data['bin_tau'])
      	  mass = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
      	  test = ((-1)*tau*np.exp((-1)*age/tau)*(tau+age)+tau*np.exp(-0/tau)*(tau+0))
      	  tempIdx = np.where(test == 0)[0]
      	  if len(tempIdx) > 0:
            t = tau[tempIdx]
            a = age[tempIdx]
      	  scale = mass/((-1)*tau*np.exp((-1)*age/tau)*(tau+age)+tau*np.exp(-0/tau)*(tau+0))
      	  sfr = scale*age*np.exp(-age/tau)
      	  temp_arr = np.log10(sfr)
        else:     
          temp_arr = np.log10(results['STAR FORMATION'].data['SFR'])

        mean, median, interval = get1DInterval(temp_arr, prob, levels)
        SFR_med.append(median)
        SFR_err.append(interval[0])

        if 'A1500' in results['DUST ATTENUATION'].data.dtype.names:
          temp_arr = results['DUST ATTENUATION'].data['A1500']
        else:
          temp_arr = results['DUST ATTENUATION'].data['A_1500']
               
        mean, median, interval = get1DInterval(temp_arr, prob, levels)
        AUV_med.append(median)
        AUV_err.append(interval[0])

        temp_arr = results['GALAXY PROPERTIES'].data['redshift']
               
        mean, median, interval = get1DInterval(temp_arr, prob, levels)
        redshift_med.append(median)
        redshift_err.append(interval[0])

        temp_arr = np.log10(results['STAR FORMATION'].data['max_stellar_age'])
               
        mean, median, interval = get1DInterval(temp_arr, prob, levels)
        age_med.append(median)
        age_err.append(interval[0])

    SFR_err=np.array(SFR_err)
    mass_err = np.array(mass_err)
    mass_med = np.array(mass_med)
    SFR_med = np.array(SFR_med)
    AUV_med = np.array(AUV_med)
    AUV_err = np.array(AUV_err)
    redshift_med = np.array(redshift_med)
    redshift_err = np.array(redshift_err)
    age_med = np.array(age_med)
    age_err = np.array(age_err)


    for i in range(len(mass_med)):
      mass_err[i] = abs(mass_med[i] - mass_err[i])
      SFR_err[i] = abs(SFR_err[i] - SFR_err[i])
      AUV_err[i] = abs(AUV_med[i] - AUV_err[i])
      redshift_err[i] = abs(redshift_med[i] - redshift_err[i])
      age_err[i] = abs(age_med[i] - age_err[i])

    SFR_err=np.transpose(SFR_err)
    mass_err = np.transpose(mass_err)
    AUV_err = np.transpose(AUV_err)
    redshift_err = np.transpose(redshift_err)
    age_err = np.transpose(age_err)

    outputDict = {'idArr': idArr, 'SFR_med':SFR_med, 'SFR_err':SFR_err, 'mass_med':mass_med, 'mass_err':mass_err, 'AUV_med':AUV_med, 'AUV_err':AUV_err, 'redshift_med':redshift_med, 'redshift_err':redshift_err, 'age_med':age_med, 'age_err':age_err}
    pickle.dump(outputDict,open(args.resultsDir+"medianDict.p","w"))
  '''
#
#  if args.salmon:
#
#
#    max_mass = np.max(mass_med)+0.01
#    if args.massLim is not None:
#      min_mass = args.massLim
#    else:
#      min_mass = np.min(mass_med)
#    #min_mass = 8.5
#    print 'np.min(mass_med): ', np.min(mass_med)
#    deltaMass = (max_mass-min_mass)/5
#    massBins = np.fromiter((x for x in np.arange(min_mass,max_mass,deltaMass)),np.float)
#    idxArr = np.zeros(len(mass_med))-1
#
#
#    mag_1500 = np.zeros(len(mass_med))
#    SFR_salmon = np.zeros(len(mass_med))
#    #find filter closest in wavelength to rest-frame 1500
#    for i in range(len(mass_med)):
#      print "inputCat[1].data['ID']", inputCat[1].data['ID']
#      print idArr[i]
#      if args.objectID:
#        tempIdx = np.where(inputCat[1].data['object_ID'] == idArr[i])[0]
#      else:
#        tempIdx = np.where(inputCat[1].data['ID'] == idArr[i])[0]
#      print 'tempIdx: ', tempIdx 
#      if len(tempIdx) > 1:
#        tempIdx = tempIdx[0]
#      print 'tempIdx: ', tempIdx 
#      lambda_1500 = 1500*(redshift_med[i]+1)
#      for f in range(len(filterInfo[1].data['filters'])):
#        if lambda_1500  > filterInfo[1].data['lambda_low'][f] and lambda_1500 <= filterInfo[1].data['lambda_high'][f]:
#          #this assumes units in nanoJy
#	  print (-2.5)*np.log10(inputCat[1].data[filterInfo[1].data['filters'][f]][tempIdx]*1E-23*1E-9)-48.58
#          mag_1500[i] = (-2.5)*np.log10(inputCat[1].data[filterInfo[1].data['filters'][f]][tempIdx]*1E-23*1E-9)-48.58
#          
#      print 'mag_1500[i]: ', mag_1500[i], redshift_med[i], age_med[i], AUV_med[i]
#      SFR_salmon[i] = np.log10(SFR_conv.Reddy_UV_to_SFR(mag_1500[i], redshift_med[i], age_med[i], auv = AUV_med[i], chab=True))
#      print 'SFR_salmon[i]: ', SFR_salmon[i]
#      print mass_med[i]
#      massBinIdx = np.where((mass_med[i] >= massBins) & (mass_med[i] < massBins+deltaMass))[0]
#      print 'massBinIdx: ', massBinIdx, len(massBinIdx)
#      if len(massBinIdx) ==1:
#        idxArr[i] = massBinIdx
#      else:
#        idxArr[i] = -1
#
#    #in bins of mass, calculate SFR and scatter
#    SFR_bins = np.zeros(len(massBins))
#    SFR_bin_err = np.zeros(len(massBins))
#    for i in range(len(massBins)):
#      print 'mass bin: ', massBins[i], np.where(idxArr == i)
#      tempIdxBin = np.where((np.isfinite(SFR_salmon) == True) & (idxArr == i))[0]
#      if len(tempIdxBin) > 0:
#     	 SFR_bins[i] = np.median(SFR_salmon[tempIdxBin])
#         if len(tempIdxBin) == 1:
#           print 'here'
#           SFR_bin_err[i] = 1E-3
#         else:
#     	   SFR_bin_err[i] = np.std(SFR_salmon[tempIdxBin])
#      else:
#         SFR_bins[i] = 0
#         SFR_bin_err[i] = 0
#
#    print 'massBins: ', massBins
#    print 'sfr_bins: ', SFR_bins
#    print 'SFR_bin_err: ', SFR_bin_err
#    tempIdx = np.where(SFR_bins > 0)[0]
#    results=poly.polyfit(massBins[tempIdx]+deltaMass/2., SFR_bins[tempIdx], 1, w=1/SFR_bin_err[tempIdx])
#    if not args.noPlots:
#      pylab.figure()
#      pylab.errorbar( massBins+deltaMass/2., SFR_bins, yerr=SFR_bin_err)
#      pylab.scatter(mass_med, SFR_salmon, linewidth=0, alpha=0.5)
#      pylab.scatter(mass_med, SFR_med, c='r', linewidth=0, alpha=0.5)
#      pylab.text(7, 3, 'y = '+str(results[1])+"x+"+str(results[0])+" scatter = "+str(np.median(SFR_bin_err)))
#      pylab.plot(massBins+deltaMass/2., poly.polyval(massBins+deltaMass/2.,results), 'k--')
#      pylab.savefig(args.resultsDir+"salmon_measurements.pdf")
#
#    #results_s = results
#    tempIdx = np.where((np.isfinite(SFR_salmon) == True) & (mass_med > min_mass))[0]
#    results2 = poly.polyfit(mass_med[tempIdx], SFR_salmon[tempIdx], 1)
#    print 'salmon fit: ', results2
#    print np.min(mass_med)
#    residuals = SFR_salmon[tempIdx]-(mass_med[tempIdx]*results2[1]+results2[0])
#    results_s = copy.deepcopy(results2)
#    scatter_salmon = np.std(residuals)
#
#    if not args.noPlots:
#      pylab.figure()
#      pylab.errorbar(SFR_salmon, SFR_med, yerr=SFR_err, fmt="o")
#      pylab.plot(SFR_salmon,SFR_salmon)
#      pylab.savefig(args.resultsDir+"SFR_comparison.pdf")
#    print 'scatter from salmon measurement: ', np.median(SFR_bin_err[np.where(SFR_bin_err > 0.001)]), scatter_salmon

  if args.santini:

    summCat = args.resultsDir+"pyp-beagle/data/BEAGLE_summary_catalogue.fits"
    data_fits = fits.open(summCat)
    idArr = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    redshift_med = np.asarray(data_fits['GALAXY PROPERTIES'].data['redshift_median'], dtype=np.float64)

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
        
    
    plt.hist(redshift_med, bins=30)
    plt.show()
    
    for i in range(len(mass_med)):
        tempIdx = (np.where(inputCat[1].data['ID'] == idArr[i]))[0][0]
#        print(tempIdx)
            
        lambda_1280 = 1200*(redshift_med[i]+1)
        lambda_2600 = 3000*(redshift_med[i]+1)
        
        temp_fnu_lester = []
        appMagArr = []
        appMagErr = []
        lambdaArr = []
        
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
#            SFR_err_santini.append(temp_sfr_err)
            Meurer_correction_arr.append(meurer_correction)

    SFR_santini = np.array(SFR_santini)
    SFR_err_santini = np.array(SFR_err_santini)

    print(SFR_santini)
    print(str(len(SFR_santini))+'/'+str(len(idArr)))


    test = np.isin(idArr, SFR_santini_ID)
    print(len(mass_med))
    print(len(mass_med[test]))  
    
    print(idArr)
    print(SFR_santini_ID)
    
    plt.scatter(mass_med[test], SFR_santini)
    plt.show()
    
    
    print(np.polyfit(mass_med[test], SFR_santini, 1))
    
        
#    print "**", SFR_santini.shape, SFR_err_santini.shape
#
#    if args.massLim is not None:
#      tempIdx = np.where((np.isfinite(SFR_santini) == True) & (mass_med > args.massLim))[0]
#    else: 
#      tempIdx = np.fromiter((x for x in range(len(mass_med))),np.int)
#    results2 = poly.polyfit(mass_med[tempIdx], SFR_santini[tempIdx], 1)#, w = 1./SFR_err_santini[tempIdx])
#    residuals = SFR_santini[tempIdx]-(mass_med[tempIdx]*results2[1]+results2[0])
#    #print 'tempIdx: ', tempIdx, np.max(np.abs(residuals)), np.std(residuals), np.max(SFR_santini)
#    ##perform 2sigma clipping 
#    #tempIdx2 = tempIdx[np.where((np.isfinite(SFR_santini[tempIdx]) == True) & (np.abs(residuals) < 2*np.std(residuals)))[0]]
#    #print 'tempIdx2: ', tempIdx2
#    #results2 = poly.polyfit(mass_med[tempIdx2], SFR_santini[tempIdx2], 1)#, w = 1./SFR_err_santini[tempIdx2])
#    #residuals = SFR_santini[tempIdx2]-(mass_med[tempIdx2]*results2[1]+results2[0])
#    print "**", len(tempIdx), len(residuals)
#    scatter_santini = np.std(residuals)
#    print scatter_santini
#    #if 1==1:
##    if not args.noPlots:
##      pylab.figure()
##      pylab.errorbar(mass_med, SFR_santini, yerr=SFR_err_santini,fmt="o")
##      pylab.scatter(mass_med, SFR_med, c='r')
##      #pylab.scatter(SFR_santini, SFR_med)
##      #pylab.scatter(AUV_med, Meurer_correction_arr)
##      #pylab.savefig(args.resultsDir+"santini_measurements_AUV.pdf")
##      pylab.savefig(args.resultsDir+"santini_measurements.pdf")
#    #sys.exit()
#    results_santini = results
#    print 'scatter from santini measurement: ', np.std(residuals)
#    
#  if args.massLim is not None:
#    massLim = args.massLim
#  else:
#    massLim = 0.
#  tempIdx = np.where(mass_med > massLim)[0]
#  results = poly.polyfit(mass_med[tempIdx], SFR_med[tempIdx], 1)
#  residuals = SFR_med[tempIdx] - poly.polyval(mass_med[tempIdx],results)
#  if not args.noPlots:
#    pylab.figure()
#    pylab.errorbar(mass_med, SFR_med, xerr=mass_err, yerr=SFR_err, fmt="o")
#    pylab.text(7, 3, 'y = '+str(results[1])+"x+"+str(results[0])+" scatter = "+str(np.std(residuals)))
#    pylab.plot(mass_med, poly.polyval(mass_med,results))
#    pylab.savefig(args.resultsDir+"medians_fitted_results.pdf")

#  print 'scatter from medians: ', np.std(residuals)
#  print results
 
#  #for salmon, the results from fitting to bin values are given while the fitted results to individual objects is given
#  #for santini
#  if args.salmon and args.santini: 
#    outputFile = args.resultsDir+"medians_results_"
#    if args.massLim is not None:
#      tempStr = str(np.round(args.massLim,2))
#      massLimStr = tempStr.replace(".","p")
#      outputFile = outputFile+massLimStr
#    f = open(outputFile+".dat","w")
#    f.write("#med_slope med_intercept med_scatter salmon_slope salmon_intercept salmon_scatter salmon_scatter_residuals "+\
#            "santini_slope santini_intercept santini_scatter_residuals \n")
#    f.write(str(results[1])+" "+str(results[0])+" "+str(np.std(residuals))+" "+\
#            str(results_s[1])+" "+str(results_s[0])+" "+str(np.median(SFR_bin_err[np.where(SFR_bin_err > 0.01)]))+" "+str(scatter_salmon)+" "+\
#            str(results2[1])+" "+str(results2[0])+" "+str(scatter_santini)+"\n")
#            #str(results_santini[1])+" "+str(results_santini[0])+" "+str(np.median(SFR_bin_err_sant))+" "+str(scatter_santini)+"\n")
#    f.close()

parser = argparse.ArgumentParser()

parser.add_argument(
    '-r', '--results-dir',
    help="Directory containing the fits.gz output files from BEAGLE fitting",
    action="store",
    type=str,
    dest="resultsDir",
    required=True
)

#parser.add_argument(
#    '--m-star',
#    help="Measure results from M_star rather than M_tot",
#    action="store_true",
#    dest="Mstar",
#    required=False,
#    default=False
#)

parser.add_argument(
    '--m-tot',
    help="Measure results from M_tot rather than M_star",
    action="store_true",
    dest="Mtot",
    required=False,
    default=False
)

parser.add_argument(
    '--delayed',
    help="Measure results from instantaneous SFR rather than SFR column in output file",
    action="store_true",
    dest="delayed",
    required=False,
    default=False
)

parser.add_argument(
  '--salmon',
  help="Measure results with Salmon+ 2015 method",
  action="store_true",
  dest="salmon",
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
  help="Input photometric catalogue, required if --salmon flag used",
  action="store",
  type=str,
  dest="inputCat",
  required=False,
  default=None
)

parser.add_argument(
  '--filter-info',
  help="file containing filter wavelengths corresponding to filters in input catalogue",
  action="store",
  type=str,
  dest="filterInfo",
  required=False,
  default=None
)

parser.add_argument(
  '--mass-lim',
  help="apply mass limit",
  action="store",
  type=float,
  dest="massLim",
  required=False,
  default=None
)

parser.add_argument(
  '--no-plots',
  help="flag to set if not wanting to output plots",
  action="store_true",
  default=False,
  dest="noPlots",
  required=False
)

parser.add_argument(
  '--object-id',
  help="the NIRCam scenarios were run using object_ID from the photometric catalogue and the CANDELS ones with ID",
  action="store_true",
  default=False,
  dest="objectID",
  required=False
)

args = parser.parse_args()


run_main(args)



'''
appMagArr, appMagErr, lambdaArr, redshift_med = [30.0, 30.1], [0.2, 0.2], [4400.0, 6000.0], 1.7
temp_sfr, meurer_correction = Santini_SFR(appMagArr, appMagErr, lambdaArr, redshift_med)
print(temp_sfr, meurer_correction)
 

fnu=5.66 # microJy
appMag1 = ((-2.5)*np.log10(fnu*1E-23*1E-9)-48.58)

fnuJy = fnu*1e-6 #Jy
appMag2 =-2.5 * np.log10(fnu*1e-6) + 8.90


appMag3 = ((-2.5)*np.log10(fnuJy*1E-23)-48.58)



print(fnu, fnuJy)
print(appMag1, appMag2, appMag3)


'''

















