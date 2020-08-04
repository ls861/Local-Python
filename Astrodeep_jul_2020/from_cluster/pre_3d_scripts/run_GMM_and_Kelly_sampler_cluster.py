import numpy as np
import matplotlib.pyplot as plt
import numpy as np

#import linmix_GMM as linmixGMM

#import linmix_truncnorm as linmix

#lester commenting out all 4 below
#import linmix_emma as linmixGMM
#import linmix as linmix
#import linmix_MH_GMM as linmix_MH_GMM
#import linmix_MH as linmix_MH
import linmix_MH_GMM_LS as linmixGMM

import warnings
warnings.filterwarnings("ignore")
#import linmix as linmix
#import pygmmis
from sklearn import mixture
import argparse
import os
import pylab
import pickle
from astropy.io import fits
import sys
import copy

parser = argparse.ArgumentParser()

parser.add_argument(
                    '-f', '--catFile',
                    help="name of folder containing BEAGLE results",
                    action="store",
                    type=str,
                    dest="folder",
                    required=True
                    )

parser.add_argument(
                    '-g', '--fit-gmm',
                    help="fit GMM to each of the objects",
                    action="store_true",
                    default=False,
                    dest="GMM",
                    required=False
                    )

parser.add_argument(
                    '-p', '--plot-gmm',
                    help="make plots of GMM fits to posteriors",
                    action="store_true",
                    default=False,
                    dest="GMMplots",
                    required=False
                    )

parser.add_argument(
                    '-r', '--run-gibbs',
                    help="run the Gibbs sampler",
                    action="store_true",
                    default=False,
                    dest="Gibbs",
                    required=False
                    )

parser.add_argument(
                    '--n-K',
                    help="how many Gaussians to use in the GMM of the x-distribution",
                    action="store",
                    type=int,
                    default=3,
                    dest="nK",
                    required=False
                    )



parser.add_argument(
                    '--n-chains',
                    help="how many chains to set running",
                    action="store",
                    type=int,
                    default=4,
                    dest="nChains",
                    required=False
                    )



parser.add_argument(
                    '--max-iter',
                    help="maximum number of iterations",
                    action="store",
                    type=int,
                    default=10000,
                    dest="maxIter",
                    required=False
                    )

parser.add_argument(
                    '--min-iter',
                    help="minimum number of iterations",
                    action="store",
                    type=int,
                    default=5000,
                    dest="minIter",
                    required=False
                    )


parser.add_argument(
                    '--n-burn',
                    help="number of burn-in iterations",
                    action="store",
                    type=int,
                    default=1000,
                    dest="nBurn",
                    required=False
                    )

parser.add_argument(
                    '--n-gauss',
                    help="how many Gaussians to fit to PDFs",
                    action="store",
                    type=int,
                    dest="nGauss",
                    required=True
                    )

#parser.add_argument(
#                    '--delayed',
#                    help="fitting to delayed SFH",
#                    action="store_true",
#                    default=False,
#                    dest="delayed",
#                    required=False
#                    )

parser.add_argument(
                    '--re-run',
                    help="if true then re-run fitting",
                    action="store_true",
                    default=False,
                    dest="reRun",
                    required=False
                    )

parser.add_argument(
                    '--sigma',
                    help="fit mass-dependence of sigma",
                    action="store_true",
                    default=False,
                    dest="sigma",
                    required=False
                    )



parser.add_argument(
                    '--subsample',
                    help="for testing - only fit the Gibbs sampler to a subset of the objects",
                    action="store",
                    type=int,
                    dest="subsample",
                    required=False,
                    default=-1
                    )

parser.add_argument(
                    '--mass-lower',
                    help="for testing - only fit to objects with mean mass higher than this value - replaces subsamp keyword if that's set, currently only works with nGauss = 1",
                    action="store",
                    type=float,
                    dest="massLowerLim",
                    required=False,
                    default=-1
                    )

parser.add_argument(
                    '--alpha-start',
                    help="starting guess for alpha",
                    action="store",
                    type=float,
                    dest="alphaStart",
                    required=False,
                    default=None
                    )

parser.add_argument(
                    '--beta-start',
                    help="starting guess for beta",
                    action="store",
                    type=float,
                    dest="betaStart",
                    required=False,
                    default=None
                    )


parser.add_argument(
                    '--no-limits',
                    help="IF --n-gauss = 1 then you can choose to not impose limits in the GMM modelling by setting this flag",
                    action="store_false",
                    dest="imposeLimits",
                    required=False,
                    default=True
                    )


parser.add_argument(
                    '--m-tot',
                    help="fit GMMs using Mtot rather than Mstar",
                    action="store_true",
                    dest="mTot",
                    required=False,
                    default=False
                    )


parser.add_argument(
                    '--delayed',
                    help="fit GMMs using instantaneous SFR",
                    action="store_true",
                    dest="delayed",
                    required=False,
                    default=False
                    )




parser.add_argument(
                    '--no-plots',
                    help="set flag if not plotting results here",
                    action="store_true",
                    dest="noPlots",
                    required=False,
                    default=False
                    )


args = parser.parse_args()

def imposeLimits(data):
    return (data[:,1] > data[:,0]-np.log10(1.185E9)) & (data[:,1] < data[:,0]-np.log10(1.E7))

fileStr = str(args.nGauss)
if not args.imposeLimits:
  fileStr = fileStr+"_noLimits"

if args.delayed:
  fileStr = fileStr+"_instSfr"

error_setup = np.geterr()
np.seterr(all="raise")

os.system("mkdir "+args.folder+"sklearn_GMM_"+fileStr+"/")

if args.GMM is True:
    x = []
    y = []
    xycov = []
    xsig = []
    ysig = []
    amp = []
    if args.reRun:
        os.system('rm '+args.folder+"sklearn_GMM_"+fileStr+"/*.p")

    fileList = os.listdir(args.folder)
    for file in fileList:
        print 'file: ', file
        if '.fits.gz' in file:
            #check if the results for this individual object already stored
            pickleName = file.replace(".fits.gz","_GMM.p")
            fileTest = os.path.isfile(args.folder+"sklearn_GMM_"+fileStr+"/"+pickleName)
            if args.reRun or fileTest == False:
                print file
                data = fits.open(args.folder+file)
                tempIdx = np.where(data['STAR FORMATION'].data['sfr'] > 0)[0]
                sfr = copy.deepcopy(data['STAR FORMATION'].data['sfr'])
                sfr[tempIdx] = np.log10(data['STAR FORMATION'].data['sfr'][tempIdx])
                tempIdx = np.where(data['STAR FORMATION'].data['sfr'] <= 0)[0]
                sfr[tempIdx] = -30
                probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
                probs_prop = probs_prop/probs_prop.sum().astype(np.float64)
                if args.delayed:
                    tau = np.power(10,data['POSTERIOR PDF'].data['tau'])
                    age = np.power(10,data['POSTERIOR PDF'].data['max_stellar_age'])
                    norm = np.power(10,data['POSTERIOR PDF'].data['mass'])/tau*np.exp((-1)*age/tau)*(tau+age)+np.power(tau,2)
                    sfr = np.log10(norm*age*np.exp((-1)*age/tau))
                if args.mTot:
                    print "definitely fitting to Mtot"
                    mass = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
                else:
                    mass = np.float64(data['GALAXY PROPERTIES'].data['M_star'])
                mass = np.log10(mass)
    #            if args.delayed:
    #                age = np.float64(data['STAR FORMATION'].data['max_stellar_age'])
    #                tau = np.float64(data['STAR FORMATION BINS'].data['bin_tau'])
    #                mTot = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
    #                print 'tau, mass, age, scale: ', np.min(tau), np.min(mTot), np.min(age)
    #                test = ((-1)*tau*np.exp((-1)*age/tau)*(tau+age)+tau*np.exp(-0/tau)*(tau+0))
    #                tempIdx = np.where(test == 0)[0]
    #                if len(tempIdx) > 0:
    #                    t = tau[tempIdx]
    #                    a = age[tempIdx]
    #                    print '**', test[tempIdx], tau[tempIdx], mTot[tempIdx], age[tempIdx], t*np.exp(-a/t)*(t+a), t*np.exp(-0/t)*(t+0)
    #                scale = mTot/((-1)*tau*np.exp((-1)*age/tau)*(tau+age)+tau*np.exp(-0/tau)*(tau+0))
    #                sfr = scale*age*np.exp(-age/tau)
    #                sfr = np.log10(sfr)
                idx = np.random.choice(len(probs_prop),size=1000, p=probs_prop)
                #gmm = pygmmis.GMM(K=1, D=2)      # K components, D dimensions
                #gmm = pygmmis.GMM(K=args.nGauss, D=2)      # K components, D dimensions
                dataIn = np.array([mass[idx],sfr[idx]])

                tryCounter = 0
                fitted = False
                while fitted == False and tryCounter <= 3:
                    tryCounter = tryCounter+1
                    #logL, U = pygmmis.fit(gmm, np.transpose(dataIn), init_method='minmax')
                    #print 'logL, U: ', logL, U
                    try:
                        if args.nGauss == 1:
                            if args.imposeLimits:
                                gmm = mixture.GaussianMixture(n_components=args.nGauss, covariance_type="full").fit(np.transpose(dataIn))
                                #logL, U = pygmmis.fit(gmm, np.transpose(dataIn), init_method='minmax', sel_callback=imposeLimits)
                            else:
                                gmm = mixture.GaussianMixture(n_components=args.nGauss, covariance_type="full").fit(np.transpose(dataIn))
                                #logL, U = pygmmis.fit(gmm, np.transpose(dataIn), init_method='minmax')
                        else:
                            gmm = mixture.GaussianMixture(n_components=args.nGauss, covariance_type="full").fit(np.transpose(dataIn))
                            #logL, U = pygmmis.fit(gmm, np.transpose(dataIn), init_method='minmax')

                        if args.nGauss > 1:
                            #Have to check that none of these give singular matrices!
                            minDet = 1E10
                            for j in range(args.nGauss):
                                cov = gmm.covariances_[j]
                                #print cov
                                #sys.exit
                                det = cov[0,0]*cov[1,1]-cov[0,1]*cov[1,0]
                                if det < minDet:
                                    minDet = det

                            if minDet > 1E-15:
                                fitted=True
                        else:
                            fitted=True
                    except:
                        print('failed ', tryCounter)


                if args.GMMplots is True:
                    obs_size = 100
                    samples,dummy = gmm.sample(n_samples = obs_size)
                    pylab.figure()
                    pylab.scatter(mass[idx],sfr[idx])
                    pylab.scatter(samples[:,0],samples[:,1],c='r')
                    plotName = file.replace(".fits.gz","_GMM_fit.pdf")
                    pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/"+plotName)
                    pylab.close()

                if fitted == True:
                  cov = gmm.covariances_
                  if args.nGauss == 1:
      #                x.append(gmm.mean[0,0])
      #                y.append(gmm.mean[0,1])#
      #                xsig.append(np.sqrt(cov[0][0][0]))
      #                ysig.append(np.sqrt(cov[0][1][1]))
      #                xycov.append(cov[0][0][1])
                      pickle.dump({'x':gmm.means_[0,0],'y':gmm.means_[0,1],'xsig':np.sqrt(cov[0][0][0]),'ysig':np.sqrt(cov[0][1][1]),'xycov':cov[0][0][1],'amp':gmm.weights_}, open(args.folder+"sklearn_GMM_"+fileStr+"/"+pickleName,'w'))
                  else:
      #                x.append(gmm.mean[:,0])
      #                y.append(gmm.mean[:,1])
                      temp_xsig = []
                      temp_ysig = []
                      temp_xycov = []
                      for j in range(args.nGauss):
                          temp_xsig.append(np.sqrt(cov[j][0][0]))
                          temp_ysig.append(np.sqrt(cov[j][1][1]))
                          temp_xycov.append(cov[j][0][1])
    #                  xsig.append(temp_xsig)
    #                  ysig.append(temp_ysig)
    #                  xycov.append(temp_xycov)
                      pickle.dump({'x':gmm.means_[:,0],'y':gmm.means_[:,1],'xsig':temp_xsig,'ysig':temp_ysig,'xycov':temp_xycov, 'amp':gmm.weights_}, open(args.folder+"sklearn_GMM_"+fileStr+"/"+pickleName,'w'))

#            amp.append(gmm.amp)

#    x = np.array(x)
#    y = np.array(y)
#    xycov = np.array(xycov)
#    xsig = np.array(xsig)
#    ysig = np.array(ysig)
#    amp = np.array(amp)
#
#    pickle.dump({'x':x,'y':y,'xycov':xycov,'xsig':xsig,'ysig':ysig,'amp':amp},open(args.folder+"sklearn_GMM_"+str(args.nGauss)+"/GMM_fits.p",'w'))

if args.Gibbs:
    print('runGibbs LESTER')
    # =============================================================================
    # LESTER REPLACE RESULTSP
    # =============================================================================

    summCat = args.folder+"pyp-beagle/data/BEAGLE_summary_catalogue.fits"
    #        summCat = "/Users/lester/Documents/PhD/Kelly/fit_106_DE/"+"pyp-beagle/data/BEAGLE_summary_catalogue.fits"
    data_fits = fits.open(summCat)
    ids = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    data = {}
    data['mass'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_tot_median'], dtype=np.float64))}
    data['mstar'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_star_median'], dtype=np.float64))}
    print('LESTER2')
    test_lester = np.asarray(data_fits['STAR FORMATION'].data['SFR_median'], dtype=np.float64)
    print('test', test_lester)
    print('asarray', np.asarray(test_lester, dtype=np.float64))
    print('LESTER2.5')
    test_lester[test_lester==0] = 1e-40
    print(min(test_lester))
    print('LESTER2.6')

    data['sfr'] = {'value':np.log10(test_lester)}
    print('LESTER3')
    data_fits.close()
    print('LESTER4')
    #        pickle.dump([data, ids], open('/Users/lester/Documents/PhD/Kelly/fit_106_DE/'+'results.p', 'wb'))
    pickle.dump([data, ids], open(args.folder+'results.p', 'wb'))
    print('LESTER5')

    #        NOTE: results.p needs to be created from a summary catalogue which contains exactly the gz files in results folder. To make something work later on, it was required that I made a subset of my summary catalogue and used that, because I was only using a subset of gz files to save time.

    # =============================================================================
    # creating input alpha beta sigma file
    # =============================================================================

    alpha_beta_scatter = {'alpha':-6.1548276141996485, 'beta':0.7560501633562806, 'scatter':0.3}
    pickle.dump(alpha_beta_scatter, open(args.folder+'alpha_beta_scatter'+str(args.nGauss)+'.p', 'wb'))
    ###pickle.dump(alpha_beta_scatter, open('/Users/lester/Documents/PhD/Kelly/fit_106_DE/'+'alpha_beta_scatter3.p', 'wb'))


    # =============================================================================
    #
    # =============================================================================




    np.seterr(under="warn")
    print args.folder+"sklearn_GMM_"+fileStr+"/GMM_fits.p"
    test = os.path.isfile(args.folder+"sklearn_GMM_"+fileStr+"/GMM_fits.p")
    print 'test: ', test
    if test == False:
        x = []
        y = []
        xsig = []
        ysig = []
        xycov = []
        amp = []
        idArr = []
        fileList = os.listdir(args.folder+"sklearn_GMM_"+fileStr)
        for file in fileList:
            print file
            if '.p' in file:
                print 'yes'
            if '.pdf' not in file:
                print 'no'
            if '.p' in file and '.pdf' not in file and 'Gibbs_results' not in file:
                data = pickle.load(open(args.folder+"sklearn_GMM_"+fileStr+"/"+file,"r"))
#                print data.keys()
                det = np.power(data['xsig'],2)*np.power(data['ysig'],2)-np.power(data['xycov'],2)
#                print 'det: ', det, data['xsig'], data['ysig'], data['xycov']
                if np.min(det) > 1E-15: #don't include any singular matrices!
                    print file
                    idArr.append(file.replace("_BEAGLE_GMM.p",""))
                    print idArr[-1]
                    x.append(data['x'])
                    y.append(data['y'])
                    xsig.append(data['xsig'])
                    ysig.append(data['ysig'])
                    xycov.append(data['xycov'])
                    amp.append(data['amp'])
        x = np.array(x)
        y = np.array(y)
        xycov = np.array(xycov)
        xsig = np.array(xsig)
        ysig = np.array(ysig)
        amp = np.array(amp)
        idArr = np.array(idArr)
        pickle.dump({'id':idArr, 'x':x,'y':y,'xycov':xycov,'xsig':xsig,'ysig':ysig,'amp':amp},open(args.folder+"sklearn_GMM_"+fileStr+"/GMM_fits.p",'w'))
        print 'dumped: ', args.folder
#    GMM = pickle.load(open('/Users/lester/Documents/PhD/Kelly/fit_106_DE/'+"sklearn_GMM_"+'3_noLimits'+"/GMM_fits.p","r"))
    GMM = pickle.load(open(args.folder+"sklearn_GMM_"+fileStr+"/GMM_fits.p","r"))
#    print GMM.keys()
#    print "GMM['xsig']: ", GMM['xsig']
#    print "GMM['ysig']: ", GMM['ysig']
#    print "GMM['xycov']: ", GMM['xycov']
    #print ' corr: ', GMM['xycov']/(GMM['xsig']*GMM['ysig'])
    det =np.power(GMM['xsig'],2)*np.power(GMM['ysig'],2)-np.power(GMM['xycov'],2)
    ####print 'det: ', np.min(det), det # LESTER this causes carnage
    #sys.exit()
#    print GMM['x'].shape
#    print 'x: ', GMM['x']
#    print 'y: ', GMM['y']
#    print 'xycov: ', GMM['xycov']
    tempIdx= np.fromiter((x for x in range(len(GMM['x']))),np.int)
    if args.subsample > 0 and args.subsample < len(GMM['x']):
        tempIdx = np.random.choice(len(GMM['x']), size=args.subsample,replace=False)
    if args.massLowerLim > 0:
        print args.folder

        data,ids = pickle.load(open(args.folder+"results.p",'r'))
        #print np.max(ids)
        tempMass = np.zeros(len(GMM['id']))

        for i in range(len(GMM['id'])):
            #print i
            #print np.max(ids)
            #print np.max(np.array(GMM['id'],np.int))

            tempIdx = np.where(ids == np.int(GMM['id'][i]))[0] # LESTER

            #print 'len(tempIdx): ', len(tempIdx)
            if len(tempIdx) == 1:
              if args.mTot:
#                print 'mTot'
                tempMass[i] = np.array(data['mass']['value'])[tempIdx]
              else:
                print 'mStar'
                tempMass[i] = np.array(data['mstar']['value'])[tempIdx]
             # print GMM['x'][i], tempMass[i]
            else:
              tempMass[i] = -1
        tempIdx = np.where(tempMass > args.massLowerLim)[0]
    print args.massLowerLim, 'len(tempIdx): ', len(tempIdx)
#        if args.nGauss == 1:
#            tempIdx = np.where(GMM['x'] > args.massLowerLim)[0]
#        else:
#            tempX = np.zeros(len(GMM['x']))
#            print len(tempX)
#            for i in range(len(tempX)):
#                tempX[i] = np.mean(GMM['x'][i])
#            tempIdx = np.where(tempX > args.massLowerLim)[0]
#            print len(tempIdx)


    if len(tempIdx) <=5:
        sys.exit()
    if args.sigma:
        plotSuffix = "_sigma"
        if args.nGauss == 1:
            print('lm = linmix_MH')
            lm = linmix_MH.LinMix(GMM['x'][tempIdx], GMM['y'][tempIdx], GMM['xsig'][tempIdx], GMM['ysig'][tempIdx], \
                                  xycov = GMM['xycov'][tempIdx], K=args.nK,nchains=args.nChains, alpha_start = args.alphaStart,\
                                  beta_start = args.betaStart)
        else:
            print('linmix_MH_GMM')
            lm = linmix_MH_GMM.LinMix(GMM['x'][tempIdx], GMM['y'][tempIdx], GMM['xsig'][tempIdx], GMM['ysig'][tempIdx], \
                                      xycovArr = GMM['xycov'][tempIdx], K=args.nK, \
                                  nGMM_err=args.nGauss, pi_err=GMM['amp'][tempIdx],nchains=args.nChains)
        lm.run_mcmc(miniter = args.minIter, maxiter = args.maxIter, n_burn=args.nBurn)
    else:
        plotSuffix = ""
        if args.nGauss == 1:
            print('lm = linmix')
            lm = linmix.LinMix(GMM['x'][tempIdx], GMM['y'][tempIdx], GMM['xsig'][tempIdx], GMM['ysig'][tempIdx], \
                               xycov = GMM['xycov'][tempIdx], K=args.nK,nchains=args.nChains)
        else:
            print "using linmixGMM"
            print('lm = linmixGMM')
            lm = linmixGMM.LinMix(GMM['x'][tempIdx], GMM['y'][tempIdx], GMM['xsig'][tempIdx], GMM['ysig'][tempIdx], \
                                  xycovArr = GMM['xycov'][tempIdx], K=args.nK, \
                                  nGMM_err=args.nGauss, pi_err=GMM['amp'][tempIdx],nchains=args.nChains)
        print 'setting the sampler running'


        massStr_lester = str(round(args.massLowerLim,2))
        massStr_lester = massStr_lester.replace(".","p")
        fold = args.folder[69]
        sbf = './linmix_inputs_{}/'.format(fold) + massStr_lester + '_' # subfolder NEEDS TO EXIST
        np.save(sbf+'id', GMM['id'][tempIdx])
        np.save(sbf+'GMMx', GMM['x'][tempIdx])
        np.save(sbf+'GMMy', GMM['y'][tempIdx])
        np.save(sbf+'GMMxsig', GMM['xsig'][tempIdx])
        np.save(sbf+'GMMysig', GMM['ysig'][tempIdx])
        np.save(sbf+'GMMxycov', GMM['xycov'][tempIdx])
        np.save(sbf+'nK', args.nK)
        np.save(sbf+'nGauss', args.nGauss)
        np.save(sbf+'pi_err', GMM['amp'][tempIdx])
        np.save(sbf+'nChains', args.nChains)
        np.save(sbf+'minIter', args.minIter)
        np.save(sbf+'maxIter', args.maxIter)

        lm.run_mcmc(miniter = args.minIter, maxiter = args.maxIter)

        np.save(sbf+'lm_chain', lm.chain)

#    lm = linmix.LinMix(GMM['x'], GMM['y'], GMM['xsig'], GMM['ysig'], xycov = GMM['xycov'], K=3)
#    lm.run_mcmc()

    print('***', len(lm.chain['alpha']))
    if args.massLowerLim > 0:
        print('LESTER TEST5') #yes
        massStr = str(round(args.massLowerLim,2))
        massStr = massStr.replace(".","p")
        plotSuffix = plotSuffix+"_massLowerLim_"+massStr


    if not args.noPlots:
      print('LESTER TEST55') #no
      pylab.figure()
      pylab.scatter(GMM['x'][tempIdx],GMM['y'][tempIdx])
      pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/chosen_object_means"+plotSuffix+".pdf")

    if args.sigma:
        print('LESTER TEST555') #no
        if not args.noPlots:
          pylab.figure()
          pylab.hist(lm.chain['sig_alpha'],bins=20)
          pylab.xlabel(r"$\sigma_{\alpha}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_alpha_hist"+plotSuffix+".pdf")
          pylab.figure()
          pylab.plot(lm.chain['sig_alpha'])
          pylab.xlabel(r"$\sigma_{\alpha}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_alpha"+plotSuffix+".pdf")
          pylab.figure()
          pylab.hist(lm.chain['sig_beta'],bins=20)
          pylab.xlabel(r"$\sigma_{\beta}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_beta_hist"+plotSuffix+".pdf")
          pylab.figure()
          pylab.plot(lm.chain['sig_beta'])
          pylab.xlabel(r"$\sigma_{\beta}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_beta"+plotSuffix+".pdf")
          pylab.figure()
          pylab.hist(lm.chain['min_sigma'],bins=20)
          pylab.xlabel(r"$\sigma_{min}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/min_sigma_hist"+plotSuffix+".pdf")
          pylab.figure()
          pylab.plot(lm.chain['min_sigma'])
          pylab.xlabel(r"$\sigma_{min}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/min_sigma"+plotSuffix+".pdf")

          pylab.figure()
          pylab.scatter(lm.chain['sig_alpha'],lm.chain['sig_beta'])
          pylab.xlabel(r"$\sigma_{\alpha}$")
          pylab.ylabel(r"$\sigma_{\beta}$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_alpha_beta"+plotSuffix+".pdf")

        outputDict = {'sig_alpha':lm.chain['sig_alpha'],\
                      'sig_beta':lm.chain['sig_beta'],\
                      'alpha':lm.chain['alpha'],\
                      'beta':lm.chain['beta']}



    else:
        print('LESTER TEST6') #yes
        if not args.noPlots:
          pylab.figure()
          pylab.hist(np.sqrt(lm.chain['sigsqr']),bins=20)
          pylab.xlabel(r"$\sigma$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma_hist"+plotSuffix+".pdf")
          pylab.figure()
          pylab.plot(np.sqrt(lm.chain['sigsqr']))
          pylab.xlabel(r"$\sigma$")
          pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/sigma"+plotSuffix+".pdf")
        outputDict = {'alpha':lm.chain['alpha'],\
                      'beta':lm.chain['beta'],\
                      #'xi':lm.chain['xi'],\
                      #'eta':lm.chain['eta'],\
                      'id':GMM['id'][tempIdx]}

        if 'xi' in lm.chain.dtype.names: #yes
          print('LESTER TEST 66')
          outputDict['xi'] = lm.chain['xi']
          outputDict['eta'] = lm.chain['eta']

    if not args.noPlots:
      print('LESTER TEST7') #no
      pylab.figure()
      pylab.hist(lm.chain['alpha'])
      pylab.xlabel(r"$\alpha$")
      pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/alpha_hist"+plotSuffix+".pdf")
      pylab.figure()
      pylab.plot(lm.chain['alpha'])
      pylab.xlabel(r"$\alpha$")
      pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/alpha"+plotSuffix+".pdf")
      pylab.figure()
      pylab.hist(lm.chain['beta'])
      pylab.xlabel(r"$\beta$")
      pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/beta_hist"+plotSuffix+".pdf")
      pylab.figure()
      pylab.plot(lm.chain['beta'])
      pylab.xlabel(r"$\beta$")
      pylab.savefig(args.folder+"sklearn_GMM_"+fileStr+"/beta"+plotSuffix+".pdf")

    print 'dumping: ', args.folder
    pickle.dump(outputDict, open(args.folder+"sklearn_GMM_"+fileStr+"/Gibbs_results"+plotSuffix+".p","w"))

np.seterr(over="warn")
