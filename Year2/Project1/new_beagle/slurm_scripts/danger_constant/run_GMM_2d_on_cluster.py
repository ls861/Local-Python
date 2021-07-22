import numpy as np
import warnings
import argparse
import os
import pickle
import sys
import copy
from sklearn import mixture
from astropy.io import fits
from scipy import integrate
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

args = parser.parse_args()

fileStr = str(args.nGauss)

if not args.mTot:
    fileStr = fileStr+"_mStar"

if args.mTot:
    fileStr = fileStr+"_mTot"

if args.delayed:
    fileStr = fileStr+"_delayed"

error_setup = np.geterr()
np.seterr(all="raise")
os.system("mkdir "+args.folder+"sklearn_GMM_2d_"+fileStr+"/")

# GMM fits
if args.reRun:
    os.system('rm '+args.folder+"sklearn_GMM_2d_"+fileStr+"/*GMM_2d.p")

fileList = os.listdir(args.folder)
for file in fileList:
    if '.fits.gz' in file:
        #check if the results for this individual object already stored
        pickleName = file.replace(".fits.gz","_GMM_2d.p")
        fileTest = os.path.isfile(args.folder+"sklearn_GMM_2d_"+fileStr+"/"+pickleName)
        if args.reRun or fileTest == False:
            print(file)
            data = fits.open(args.folder+file)

            if args.mTot:
                print "definitely fitting to Mtot"
                mass = np.float64(data['GALAXY PROPERTIES'].data['M_tot'])
            else:
                print "definitely fitting to Mstar"
                mass = np.float64(data['GALAXY PROPERTIES'].data['M_star'])
            mass = np.log10(mass)

            if not args.delayed:
                temp_sfr = np.float64(data['STAR FORMATION'].data['SFR'])
                sfr = np.log10(np.where(temp_sfr<=0, 1e-30, temp_sfr))
            else:
                ID = file.replace('_BEAGLE.fits.gz','')
                sfr = np.load(args.folder+'{}_log_instant_sfr.npy'.format(ID))

            probs_prop = np.array(data['POSTERIOR PDF'].data['probability'], np.float64)
            probs_prop = probs_prop/probs_prop.sum().astype(np.float64)

            redshift = np.float64(data['GALAXY PROPERTIES'].data['redshift'])

            idx = np.random.choice(len(probs_prop), size=1000, p=probs_prop)
            dataIn = np.array([mass[idx],sfr[idx]])

            tryCounter = 0
            fitted = False
            while fitted == False and tryCounter <= 50:
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
                temp_xycov = []

                for j in range(args.nGauss):
                    temp_xsig.append(np.sqrt(cov[j][0][0]))
                    temp_ysig.append(np.sqrt(cov[j][1][1]))
                    temp_xycov.append(cov[j][0][1])

                pickle.dump({'x':gmm.means_[:,0], 'y':gmm.means_[:,1], 'xsig':temp_xsig, 'ysig':temp_ysig, 'xycov':temp_xycov, 'amp':gmm.weights_}, open(args.folder+"sklearn_GMM_2d_"+fileStr+"/"+pickleName,'w'))

# Create Kelly input from GMM fits

fold = args.folder[85]
if not args.mTot:
    fold = fold+"_mStar"
if args.mTot:
    fold = fold+"_mTot"
if args.delayed:
    fold = fold+"_delayed"

os.system('mkdir ./kelly_2d_GMM_inputs/')
sbf = './kelly_2d_GMM_inputs/{}_'.format(fold)

if args.reRun:
    os.system('rm '+sbf+'GMM_2d_fits.p')

np.seterr(under="warn")
fileTest = os.path.isfile(sbf+'GMM_2d_fits.p')

if args.reRun or fileTest == False:

    idArr = []
    x = []
    y = []
    xsig = []
    ysig = []
    xycov = []
    amp = []

    fileList = os.listdir(args.folder+"sklearn_GMM_2d_"+fileStr)

    for file in fileList:
        print(file)
        if '_2d.p' in file and '.pdf' not in file:
            data = pickle.load(open(args.folder+"sklearn_GMM_2d_"+fileStr+"/"+file,"r"))

            # lester
            idArr.append(file.replace("_BEAGLE_GMM_2d.p",""))
            x.append(data['x'])
            y.append(data['y'])
            xsig.append(data['xsig'])
            ysig.append(data['ysig'])
            xycov.append(data['xycov'])
            amp.append(data['amp'])

    idx_sort = np.argsort(np.asarray(np.array(idArr), float))

    idArr = np.array(idArr)[idx_sort]
    x = np.array(x)[idx_sort]
    y = np.array(y)[idx_sort]
    xycov = np.array(xycov)[idx_sort]
    xsig = np.array(xsig)[idx_sort]
    ysig = np.array(ysig)[idx_sort]
    amp = np.array(amp)[idx_sort]

    pickle.dump({'id':idArr, 'x':x, 'y':y, 'xycov':xycov, 'xsig':xsig, 'ysig':ysig, 'amp':amp, 'nGauss':args.nGauss},open(sbf+'GMM_2d_fits.p','w'))
    print('dumped: ', args.folder)

    os.system('rm -r '+args.folder+"sklearn_GMM_2d_"+fileStr+"/")

np.seterr(over="warn")
