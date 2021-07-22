import numpy as np
import warnings
import argparse
import os
import pickle
import sys
import copy
from sklearn import mixture
from astropy.io import fits
warnings.filterwarnings("ignore")

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
                    '--n-gauss',
                    help="how many Gaussians to fit to PDFs",
                    action="store",
                    type=int,
                    dest="nGauss",
                    required=True
                    )

parser.add_argument(
                    '--re-run',
                    help="if true then re-run fitting",
                    action="store_true",
                    default=False,
                    dest="reRun",
                    required=False
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
                    '--m-tot',
                    help="fit GMMs using Mtot rather than Mstar",
                    action="store_true",
                    dest="mTot",
                    required=False,
                    default=False
                    )

parser.add_argument(
                    '--no-limits',
                    help="IF --n-gauss = 1 then you can choose to not impose limits in the GMM modelling by setting this flag",
                    action="store_false",
                    dest="imposeLimits",
                    required=False,
                    default=True
                    )

args = parser.parse_args()

fileStr = str(args.nGauss)
if not args.imposeLimits:
  fileStr = fileStr+"_noLimits"
error_setup = np.geterr()
np.seterr(all="raise")
os.system("mkdir "+args.folder+"sklearn_GMM_"+fileStr+"/")

if args.GMM is True:
    if args.reRun:
        os.system('rm '+args.folder+"sklearn_GMM_"+fileStr+"/*.p")

    fileList = os.listdir(args.folder)
    for file in fileList:
        if '.fits.gz' in file:
            #check if the results for this individual object already stored
            pickleName = file.replace(".fits.gz","_GMM_3d.p")
            fileTest = os.path.isfile(args.folder+"sklearn_GMM_"+fileStr+"/"+pickleName)
            if args.reRun or fileTest == False:
                print(file)
                data = fits.open(args.folder+file)
                tempIdx = np.where(data['STAR FORMATION'].data['sfr'] > 0)[0]
                sfr = copy.deepcopy(data['STAR FORMATION'].data['sfr'])
                sfr[tempIdx] = np.log10(data['STAR FORMATION'].data['sfr'][tempIdx])
                tempIdx = np.where(data['STAR FORMATION'].data['sfr'] <= 0)[0]
                sfr[tempIdx] = -30
                probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
                probs_prop = probs_prop/probs_prop.sum().astype(np.float64)

                if args.mTot:
                    print "definitely fitting to Mtot"
                    mass = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
                else:
                    mass = np.float64(data['GALAXY PROPERTIES'].data['M_star'])
                mass = np.log10(mass)

                redshift = np.float64(data['GALAXY PROPERTIES'].data['redshift'])

                idx = np.random.choice(len(probs_prop),size=1000, p=probs_prop)
                dataIn = np.array([mass[idx],sfr[idx],redshift[idx]])

                tryCounter = 0
                fitted = False
                while fitted == False and tryCounter <= 10:
                    tryCounter = tryCounter+1

                    try:
                        print('fitting GMM')
                        gmm = mixture.GaussianMixture(n_components=args.nGauss, covariance_type="full").fit(np.transpose(dataIn))

                        #Have to check that none of these give singular matrices!
                        minDet = 1E10
                        for j in range(args.nGauss):
                            cov = gmm.covariances_[j]
                            det = np.linalg.det(cov)

                            if det < minDet:
                                minDet = det

                        if minDet > 1E-15:
                            fitted=True

                    except:
                        print('failed ', tryCounter)

                if fitted == True:
                    cov = gmm.covariances_

                    temp_xsig = []
                    temp_ysig = []
                    temp_zsig = []
                    temp_xycov = []
                    temp_xzcov = []
                    temp_yzcov = []

                    for j in range(args.nGauss):
                        temp_xsig.append(np.sqrt(cov[j][0][0]))
                        temp_ysig.append(np.sqrt(cov[j][1][1]))
                        temp_zsig.append(np.sqrt(cov[j][2][2]))
                        temp_xycov.append(cov[j][0][1])
                        temp_xzcov.append(cov[j][0][2])
                        temp_yzcov.append(cov[j][1][2])

                    pickle.dump({'x':gmm.means_[:,0], 'y':gmm.means_[:,1], 'z':gmm.means_[:,2], 'xsig':temp_xsig, 'ysig':temp_ysig, 'zsig':temp_zsig, 'xycov':temp_xycov, 'xzcov':temp_xzcov, 'yzcov':temp_yzcov, 'amp':gmm.weights_}, open(args.folder+"sklearn_GMM_"+fileStr+"/"+pickleName,'w'))

if args.Gibbs:

    summCat = args.folder+"pyp-beagle/data/BEAGLE_summary_catalogue.fits"
    data_fits = fits.open(summCat)
    ids = np.asarray(data_fits['POSTERIOR PDF'].data['ID'], dtype=int)
    print(ids)
    data = {}
    data['mass'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_tot_median'], dtype=np.float64))}
    data['mstar'] = {'value':np.log10(np.asarray(data_fits['GALAXY PROPERTIES'].data['M_star_median'], dtype=np.float64))}
    temp_sfr = np.asarray(data_fits['STAR FORMATION'].data['SFR_median'], dtype=np.float64)
    data['sfr'] = {'value':np.log10( np.where(temp_sfr==0, 1e-50, temp_sfr) )}
    data_fits.close()
    pickle.dump([data, ids], open(args.folder+'results.p', 'wb'))

    np.seterr(under="warn")
    test = os.path.isfile(args.folder+"sklearn_GMM_"+fileStr+"/GMM_3d_fits.p")

    if test == False:

        idArr = []
        x = []
        y = []
        z = []
        xsig = []
        ysig = []
        zsig = []
        xycov = []
        xzcov = []
        yzcov = []
        amp = []

        fileList = os.listdir(args.folder+"sklearn_GMM_"+fileStr)

        for file in fileList:
            print(file)
            if '_3d.p' in file:
                print 'yes'
            if '.pdf' not in file:
                print 'no'
            if '_3d.p' in file and '.pdf' not in file:
                data = pickle.load(open(args.folder+"sklearn_GMM_"+fileStr+"/"+file,"r"))

                # lester
                idArr.append(file.replace("_BEAGLE_GMM_3d.p",""))
                x.append(data['x'])
                y.append(data['y'])
                z.append(data['z'])
                xsig.append(data['xsig'])
                ysig.append(data['ysig'])
                zsig.append(data['zsig'])
                xycov.append(data['xycov'])
                xzcov.append(data['xzcov'])
                yzcov.append(data['yzcov'])
                amp.append(data['amp'])

        idArr = np.array(idArr)
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        xycov = np.array(xycov)
        xzcov = np.array(xzcov)
        yzcov = np.array(yzcov)
        xsig = np.array(xsig)
        ysig = np.array(ysig)
        zsig = np.array(zsig)
        amp = np.array(amp)

        pickle.dump({'id':idArr, 'x':x, 'y':y, 'z':z, 'xycov':xycov, 'xzcov':xzcov, 'yzcov':yzcov, 'xsig':xsig, 'ysig':ysig, 'zsig':zsig, 'amp':amp},open(args.folder+"sklearn_GMM_"+fileStr+"/GMM_3d_fits.p",'w'))
        print 'dumped: ', args.folder

    GMM = pickle.load(open(args.folder+"sklearn_GMM_"+fileStr+"/GMM_3d_fits.p","r"))

    tempIdx= np.fromiter((x for x in range(len(GMM['x']))),np.int)

    if args.massLowerLim > 0:

        data,ids = pickle.load(open(args.folder+"results.p",'r'))
        tempMass = np.zeros(len(GMM['id']))

        for i in range(len(GMM['id'])):

            tempIdx = np.where(ids == np.int(GMM['id'][i]))[0]

            if len(tempIdx) == 1:
              if args.mTot:
                tempMass[i] = np.array(data['mass']['value'])[tempIdx]
              else:
                print 'mStar'
                tempMass[i] = np.array(data['mstar']['value'])[tempIdx]
            else:
              tempMass[i] = -1
        tempIdx = np.where(tempMass > args.massLowerLim)[0]

    if len(tempIdx) <=5:
        sys.exit()

    plotSuffix = ''

    massStr_lester = str(round(args.massLowerLim,2))
    massStr_lester = massStr_lester.replace('.','p')

    fold = args.folder[69]
    os.system('mkdir ./linmix_inputs_GMM_3d_{}/'.format(fold))
    sbf = './linmix_inputs_GMM_3d_{}/{}_'.format(fold,massStr_lester)

    np.save(sbf+'id', GMM['id'][tempIdx])
    np.save(sbf+'GMMx', GMM['x'][tempIdx])
    np.save(sbf+'GMMy', GMM['y'][tempIdx])
    np.save(sbf+'GMMz', GMM['z'][tempIdx])
    np.save(sbf+'GMMxsig', GMM['xsig'][tempIdx])
    np.save(sbf+'GMMysig', GMM['ysig'][tempIdx])
    np.save(sbf+'GMMzsig', GMM['zsig'][tempIdx])
    np.save(sbf+'GMMxycov', GMM['xycov'][tempIdx])
    np.save(sbf+'GMMxzcov', GMM['xzcov'][tempIdx])
    np.save(sbf+'GMMyzcov', GMM['yzcov'][tempIdx])
    np.save(sbf+'nK', args.nK)
    np.save(sbf+'nGauss', args.nGauss)
    np.save(sbf+'pi_err', GMM['amp'][tempIdx])
    np.save(sbf+'nChains', args.nChains)
    np.save(sbf+'minIter', args.minIter)
    np.save(sbf+'maxIter', args.maxIter)

np.seterr(over="warn")
