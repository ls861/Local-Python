import numpy as np
from astropy.io import fits
import argparse
import pickle
from pathos.multiprocessing import ProcessingPool
import os
from astropy.cosmology import FlatLambdaCDM
import subprocess
#import multi_level_modelling_truncated_sigma_prop_halfnorm_prior as mod
import warnings
warnings.filterwarnings("ignore")
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)



def run_input_vs_output(inputDict):
  #os.system("python -W ignore input_vs_output_plots.py -r "+inputDict['folder']+" -c "+inputDict['inputCat']+" --force-sfr --force-mstar --rotated --start-id 1 --no-plots --sym-link")
  os.system("python -W ignore input_vs_output_plots.py -r "+inputDict['folder']+" -c "+inputDict['inputCat']+" --force-sfr --force-mstar --start-id "+str(inputDict['startId'])+" --no-plots --sym-link")

def run_GMM_fits(inputDict):
  command = "python -W ignore run_GMM_and_Kelly_sampler.py -f "+inputDict['folder']+" -g --n-gauss "+inputDict['nGauss']+" --no-limits"
  if inputDict['mTot']:
    command = command+" --m-tot"
  if inputDict['delayed']:
    command = command+" --delayed"
  print command
  os.system(command)


def run_Kelly07(inputDict):
  print 'here'
  success = False
  counter = 1
  while success == False and counter < 5:
    counter = counter+1
    command = "python -W ignore run_GMM_and_Kelly_sampler.py -f "+inputDict['folder']+\
          " -r --n-gauss "+inputDict['nGauss']+" --no-limits --mass-lower "+inputDict['massLim'] +\
          " --n-chains 2 --min-iter 30000 --max-iter 30000 --no-plots --n-K 3"
    if inputDict['mTot']:
      command = command + ' --m-tot'
#    command = "python -W ignore run_GMM_and_Kelly_sampler.py -f "+inputDict['folder']+\
#          " -r --n-gauss "+inputDict['nGauss']+" --no-limits --mass-lower "+inputDict['massLim'] +\
#          " --n-chains 2 --min-iter 5000 --max-iter 5000 --no-plots --n-K 3 --m-tot"
    if inputDict['delayed']:
      command = command+" --delayed"
    print command
    try:
      os.system(command)
      success = True
    except:
      print('run_Kelly07 failed', counter)

def run_Kelly07_sigma_vary(inputDict):
  print 'here'
  success = False
  counter = 1
  while success == False and counter < 5:
    counter = counter+1
    try:
      #os.system("python -W ignore run_GMM_and_Kelly_sampler.py -f "+inputDict['folder']+\
      #        " -r --n-gauss "+inputDict['nGauss']+" --no-limits --mass-lower "+inputDict['massLim'] +\
      #        " --n-chains 2 --min-iter 5000 --max-iter 5000 --sigma --no-plots")
      os.system("python -W ignore run_GMM_and_Kelly_sampler.py -f "+inputDict['folder']+\
              " -r --n-gauss "+inputDict['nGauss']+" --no-limits --mass-lower "+inputDict['massLim'] +\
              " --n-chains 2 --min-iter 30000 --max-iter 50000 --n-burn 5000 --sigma")
      success = True
    except:
      print('run_uelly07_sigma_vary failed', counter)

def run_modelling_massLims(inputDict):
  n_iter = 10000
  n_burn = 5000
  #n_iter = 200
  #n_burn = 100
  beta0_in = -8.
  beta1_in = 1.
  sigma_in = 0.3
  #careful changing the 3rd element (proposal distribution width for sigma) as it is not optimised after the iterations
  tmpArr=[0.001,0.001,0.008]
  print inputDict.keys()
  print 'inputDict["massLimit"]: ', inputDict['massLimit']
  tempStr = str(inputDict['massLimit'])
  massStr = tempStr.replace(".","p")
  fileName = inputDict['folder']+"tracetest1_mLim_"+massStr+".fits"
  fileTest = os.path.isfile(fileName)
  lim = inputDict['massLimit']
  print 'Im here: ', fileName
  mStar = True
  if inputDict['mTot'] == True:
    mStar = False
  print 'mStar: ', mStar
#Horrible Emma tweak here
  if 1==1:
#  if not fileTest:
    print 'True'
    mod.run_multi_level_modelling(inputDict['folder'], n_iter, n_burn, beta0_in, beta1_in, sigma_in, nRuns=1, plotSuffix='test1', ageFree=True, sfrFree=inputDict['sfrFree'], metallicityFree=False, dustFree=True, initialProposal=tmpArr, massLimit = round(lim,2), mStar = mStar)
    mod.run_multi_level_modelling(inputDict['folder'], n_iter, n_burn, beta0_in, beta1_in, sigma_in, nRuns=1, plotSuffix='test2', ageFree=True, sfrFree=inputDict['sfrFree'], metallicityFree=False, dustFree=True, initialProposal=tmpArr, massLimit = round(lim,2), mStar = mStar)

def simplify_fits(inputDict):
  os.system("mkdir "+inputDict['folder']+"fitsFiles/")
  command = "python -W ignore simplify_fits_files.py -r "+inputDict['folder']+" -o "+inputDict['folder']+"fitsFiles/"
  if not inputDict['sfrFree']:
    command = command + " --add-sfr"
  os.system(command)

def measure_medians(inputDict):
  command = "python -W ignore M_SFR_fit_medians.py -r "+inputDict['folder']+" --salmon --input-cat "+inputDict['inputPhotCat']+\
              " --filter-info "+inputDict['filterInfo']+" --santini --mass-lim "+inputDict['massLim']
              #" --filter-info "+inputDict['filterInfo']+" --santini --no-plots --mass-lim "+inputDict['massLim']
  if not inputDict['mTot']:
    command = command+' --m-star'
  if inputDict['objectID']:
    command = command+" --object-id"
  print command
  os.system(command)

# 31 July 2018 - Emma Curtis Lake
# Script to run all analysis procedures on a given scenario - probably there will have to be some
# preparation before running this script
# 1 - Derive medians and confidence contours from the BEAGLE fits
# 2 - Derive mass limits for different fractions of objects with constraints
#     far from the limits in the priors
# 3 - Fit GMMs to pdfs, once with 1 Gaussian, once with 6 Gaussians
# 4 - For each mass limit, measure alpha,beta,sigma with:
#     a - Kelly07 algorithm
#     b - David's algorithm
#     c - Salmon
#     d - Santini
#     e - Medians
# 5 - Measure alpha,beta,sigma_alpha,sigma_beta with modified Kelly07 algorithm


#TODO - allow for flags to turn each measurement on/off
parser = argparse.ArgumentParser()

parser.add_argument(
                    '-i', '--input-cat-file',
                    help="name of catalogue containing input values, give list of files in .txt file if more than one to be used to derive the mass limits (e.g. with the test scenarios)",
                    action="store",
                    type=str,
                    dest="inputCatFile",
                    required=True
                    )

parser.add_argument(
                    '-f', '--results-folder',
                    help="name of folder containing BEAGLE fits, give list of files in .txt file if more than one to be used to derive the mass limits (e.g. with the test scenarios)",
                    action="store",
                    type=str,
                    dest="resultsFolder",
                    required=True
                    )

parser.add_argument(
                    '-i2', '--input-cat-file-2',
                    help = "name of second input catalogue if it exists (e.g. for catalogues based on two different snapshots)",
                    action="store",
                    type=str,
                    dest="inputCatFile2",
                    default=None,
                    required=False
                    )

parser.add_argument(
                    '-f2', '--results-folder-2',
                    help = "name of second results folder if it exists",
                    action="store",
                    type=str,
                    dest="resultsFolder2",
                    default=None,
                    required=False
                    )

parser.add_argument(
                    '-fj', '--results-folder-joint',
                    help = "name of joint results folder",
                    action="store",
                    type=str,
                    dest="resultsFolderJoint",
                    default=None,
                    required=False
                    )

parser.add_argument(
                    '-z', '--max-redshift',
                    help = "maxmimum redshift within sample",
                    action="store",
                    type=float,
                    dest="zMax",
                    required=True
                    )

parser.add_argument(
                    '--m-tot',
                    help = "measurements with mTot rather than mStar",
                    action="store_true",
                    default=False,
                    dest="mTot",
                    required=False
                    )


parser.add_argument(
                    '--sfr-free',
                    help = "fitted with sfr as free parameter - set flag if true",
                    action="store_true",
                    default=False,
                    dest="sfrFree",
                    required=False
                    )

parser.add_argument(
                    '-p', '--input-phot-cat',
                    help = "input photometric catalogue",
                    action="store",
                    type=str,
                    dest="inputPhotCat",
                    required=True
                    )


parser.add_argument(
                    '-filt', '--filter-info-file',
                    help = "filter information file",
                    action="store",
                    type=str,
                    dest="filterInfo",
                    required=True
                    )


parser.add_argument(
                    '--run-all',
                    help = "set flag if to run everything",
                    action="store_true",
                    default=False,
                    dest="runAll",
                    required=False
                    )

parser.add_argument(
                    '--run-input-vs-output-prep',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runInputVsOutputPrep",
                    required=False
                    )

parser.add_argument(
                    '--run-input-vs-output',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runInputVsOutput",
                    required=False
                    )

parser.add_argument(
                    '--run-measure-mass-limits',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runMeasureMassLimits",
                    required=False
                    )

parser.add_argument(
                    '--run-GMM-fits',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runGMMfits",
                    required=False
                    )

parser.add_argument(
                    '--run-Kelly07',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runKelly07",
                    required=False
                    )

parser.add_argument(
                    '--run-simplify-fits',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runSimplifyFits",
                    required=False
                    )

parser.add_argument(
                    '--run-Davids',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runDavids",
                    required=False
                    )


parser.add_argument(
                    '--run-measure-medians',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runMeasureMedians",
                    required=False
                    )


parser.add_argument(
                    '--run-Kelly07-sigma',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="runKelly07sigma",
                    required=False
                    )

parser.add_argument(
                    '--re-run',
                    help = "set flag if to run",
                    action="store_true",
                    default=False,
                    dest="reRun",
                    required=False
                    )

parser.add_argument(
                    '--id-flip',
                    help = "if fitted ids in beagle from object_ID column in phot cat",
                    action="store_true",
                    default=False,
                    dest="idFlip",
                    required=False
                    )

parser.add_argument(
                    '--regular-mlim-grid',
                    help = "set to true if you want a regular mlim grid",
                    action="store_true",
                    default=False,
                    dest="regularMlimGrid",
                    required=False
                    )


parser.add_argument(
                    '--regular-mlim-grid-high',
                    help = "set to true if you want a regular mlim grid",
                    action="store_true",
                    default=False,
                    dest="regularMlimGridHigh",
                    required=False
                    )

parser.add_argument(
                    '--delayed',
                    help = "fit GMM and Kelly07 to instantaneous SFRs for a delayed history",
                    action="store_true",
                    default=False,
                    dest="delayed",
                    required=False
                    )

parser.add_argument(
                    '--split-gmm-z',
                    help = "split up GMM.p files according by redshift bin for running the Kelly07 sampler",
                    action="store_true",
                    default=False,
                    dest="splitGMMbyRedshift",
                    required=False
                    )

parser.add_argument(
                    '--regular-z-grid',
                    help = "split up results by grid in redshift",
                    action="store_true",
                    default=False,
                    dest="regularZgrid",
                    required=False
                    )

parser.add_argument(
                    '--start-id',
                    help = "the starting ID of the catalogue",
                    action="store",
                    default=1,
                    type=int,
                    dest="startId",
                    required=False
                    )

args = parser.parse_args()

x  = np.arange(1,1000,0.1)
y  = np.arange(1,1000,0.1)


pool = ProcessingPool(nodes=6)
#pool = ProcessingPool(nodes=2)

if args.runAll:
  args.runInputVsOutput = True
  args.runMeasureMassLimits = True
  args.runGMMfits = True
  args.runKelly07 = True
  args.runSimplifyFits = True
  args.runDavids = True
  args.runMeasureMedians = True
  args.runKelly07sigma = True


fileList = False
if '.txt' in args.resultsFolder:
  folderListData = np.genfromtxt(args.resultsFolder,dtype=None, names=True)
  catListData = np.genfromtxt(args.inputCatFile, dtype=None, names=True)
  inputPhotCatData = np.genfromtxt(args.inputPhotCat, dtype=None, names=True)
  fileList = True
else:
  if args.resultsFolder2 is None:
    args.resultsFolderJoint = args.resultsFolder
  else:
    if args.resultsFolderJoint is None:
      print "error, must give resultsFolderJoint as input"
      sys.exit()

#remove any empty fits files
if fileList:
  print folderListData['folder']
  for i in range(len(folderListData['folder'])):
    temp = subprocess.check_output("ls -lth "+folderListData['folder'][i]+"*.fits.gz | grep ' 0 ' | awk '{print $9}'", shell=True)
    print 'temp: ', temp.split(".gz")
    if len(temp) > 0:
      for file in temp.split(".gz"):
        tempStr = file.replace(".fits",".fits.gz")
        tempStr2 = tempStr.replace("\n","")
        print 'removing file', tempStr2
        os.system("rm "+tempStr2)
else:
  temp = subprocess.check_output("ls -lth "+args.resultsFolderJoint+"*.fits.gz | grep ' 0 ' | awk '{print $9}'", shell=True)
  if len(temp) > 0:
    for file in temp.split(".gz"):
      tempStr = file.replace(".fits",".fits.gz")
      tempStr2 = tempStr.replace("\n","")
      print 'removing file', tempStr2
      os.system("rm "+tempStr2)

#only if flag turned on, run input_vs_output prep (so that it's done after files have been deleted that are added to
if args.runInputVsOutputPrep:
  os.system('python input_vs_output_plot_prep.py')


#first measure the output values
if args.runInputVsOutput:
  #I don't just run pyp-beagle because I need the rotated values as well!
  if fileList:
    inputArr = []
    for i in range(len(folderListData['folder'])):
      inputArr.append({'folder':folderListData['folder'][i],'inputCat':catListData['inputCat'][i],'startId':args.startId})

    pool.map(run_input_vs_output, inputArr)
  else:
    if args.inputCatFile2 is not None: #then we're working with sym links, not great logic here but there we go
      print 'symLink'
      #os.system("python -W ignore input_vs_output_plots.py -r "+args.resultsFolder+" -c "+args.inputCatFile+" --force-sfr --force-mstar --rotated --start-id 1 --no-plots --sym-link")
    else:
      print 'not symLink'
      #os.system("python -W ignore input_vs_output_plots.py -r "+args.resultsFolder+" -c "+args.inputCatFile+" --force-sfr --force-mstar --rotated --start-id 1 --no-plots")
      print "python -W ignore input_vs_output_plots.py -r "+args.resultsFolder+" -c "+args.inputCatFile+" --force-sfr --force-mstar --start-id "+str(args.startId)+" --no-plots"
      os.system("python -W ignore input_vs_output_plots.py -r "+args.resultsFolder+" -c "+args.inputCatFile+" --force-sfr --force-mstar --start-id "+str(args.startId)+" --no-plots --force-redshift")


  #because of the way input_vs_output_plots works, it requires a single output directory for a single input catalogue.
  #I have a different input catalogue for different snapshots, but made a single catalogue for BEAGLE to fit to that
  #includes objects from two different snapshots.  Therefore there are two separate folders containing symlinks
  #to the BEAGLE results split according to input catalogue.
  if args.inputCatFile2 is not None:
    #os.system("python -W ignore input_vs_output_plots.py -r "+args.resultsFolder2+" -c "+args.inputCatFile2+" --force-sfr --force-mstar --rotated --start-id 1 --no-plots --sym-link")
    os.system("python -W ignore input_vs_output_plots.py -r "+args.resultsFolder2+" -c "+args.inputCatFile2+" --force-sfr --force-mstar --start-id "+str(args.startId)+" --no-plots --sym-link")
#    if args.idFlip:
#      print 'idFlip'
#      os.system("python collate_results.py -f "+args.resultsFolder+" -f2 "+args.resultsFolder2+" -fj "+args.resultsFolderJoint+" -c "+args.inputPhotCat+" --start-id 1")
#    else:
#      os.system("python collate_results.py -f "+args.resultsFolder+" -f2 "+args.resultsFolder2+" -fj "+args.resultsFolderJoint)

#next derive mass limits
if args.runMeasureMassLimits and not args.regularMlimGrid:
  x_sfr_m_rot = []
  x_sfr_m_rot_err = []
  y_sfr_m_rot = []
  y_sfr_m_rot_err = []
  mstar = []
  if fileList:
    for i in range(len(folderListData['folder'])):
      data,ids = pickle.load(open(folderListData['folder'][i]+"results.p",'r'))
      x_sfr_m_rot.extend(data['x_sfr_m_rot']['value'])
      x_sfr_m_rot_err.extend(data['x_sfr_m_rot']['err'])
      y_sfr_m_rot.extend(data['y_sfr_m_rot']['value'])
      y_sfr_m_rot_err.extend(data['y_sfr_m_rot']['err'])
      mstar.extend(data['mstar']['value'])
  else:
    data,ids = pickle.load(open(args.resultsFolder+"results.p",'r'))
    x_sfr_m_rot.extend(data['x_sfr_m_rot']['value'])
    x_sfr_m_rot_err.extend(data['x_sfr_m_rot']['err'])
    y_sfr_m_rot.extend(data['y_sfr_m_rot']['value'])
    y_sfr_m_rot_err.extend(data['y_sfr_m_rot']['err'])
    mstar.extend(data['mstar']['value'])
  if args.resultsFolder2 is not None:
    data,ids = pickle.load(open(args.resultsFolder2+"results.p",'r'))
    x_sfr_m_rot.extend(data['x_sfr_m_rot']['value'])
    x_sfr_m_rot_err.extend(data['x_sfr_m_rot']['err'])
    y_sfr_m_rot.extend(data['y_sfr_m_rot']['value'])
    y_sfr_m_rot_err.extend(data['y_sfr_m_rot']['err'])
    mstar.extend(data['mstar']['value'])
  x_sfr_m_rot = np.array(x_sfr_m_rot)
  x_sfr_m_rot_err = np.array(x_sfr_m_rot_err)
  y_sfr_m_rot = np.array(y_sfr_m_rot)
  y_sfr_m_rot_err = np.array(y_sfr_m_rot_err)
  mstar = np.array(mstar)

  max_age = cosmo.age(args.zMax).value-cosmo.age(15).value
  limits = [0.7071*np.log10(1.E6),0.7071*np.log10(max_age*1E9)]


  deltaM = 0.1
  Mbins = np.arange(7.5,11.5,deltaM)
  fracArr = np.zeros(len(Mbins))
  for i,M in enumerate(Mbins):
    tempIdx = np.where(mstar > M)[0]
    nInBin = np.float(len(tempIdx))
    if nInBin > 10:
      tempIdx = np.where((mstar > M) & (x_sfr_m_rot-x_sfr_m_rot_err[:,0] > limits[0]+0.3) & \
                                       (x_sfr_m_rot+x_sfr_m_rot_err[:,1] < limits[1]-0.3))[0]
      fracArr[i] = np.float(len(tempIdx))/nInBin


  print Mbins, fracArr
  maxMidx = np.where(fracArr > 0)[0]
  tempIdx = np.where(Mbins < Mbins[maxMidx[-1]])[0]
  limitArr = [0.2,0.4,0.6,0.75,0.8,0.85,0.9,0.95]
  f = np.interp(limitArr, fracArr[tempIdx], Mbins[tempIdx]+deltaM/2.)
  print limitArr, f

  massLimits = {}
  for i in range(len(limitArr)):
    tempStr = str(limitArr[i])
    limitStr = tempStr.replace(".","p")
    massLimits[limitStr] = {'frac':limitArr[i],'massLimit':f[i]}

  pickle.dump(massLimits, open(args.resultsFolderJoint+"massLimits.p","w"))

#Next perform GMM fits
if args.runGMMfits:
  if fileList:
    inputArr = []
    for i in range(len(folderListData['folder'])):
      #inputArr.append({'folder':folderListData['folder'][i], 'nGauss': str(1), 'mTot': args.mTot})
      inputArr.append({'folder':folderListData['folder'][i], 'nGauss': str(3), 'mTot': args.mTot, 'delayed':args.delayed})

    pool.map(run_GMM_fits, inputArr)
  else:
    #command = "python -W ignore run_GMM_and_Kelly_sampler.py -f "+args.resultsFolderJoint+" -g --n-gauss 1 --no-limits"
    #if args.mTot:
    #  command = command + " --m-tot"
    #print command
    #sys.exit()
    #os.system(command)
    command = "python -W ignore run_GMM_and_Kelly_sampler.py -f "+args.resultsFolderJoint+" -g --n-gauss 3 --no-limits"
    if args.mTot:
      command = command + " --m-tot"
    os.system(command)

if args.splitGMMbyRedshift:
  zArr = [1.5,2.5,3.5,4.5,5.5]
  deltaZ = 1.
  zStrArr = ['1p5','2p5','3p5','4p5','5p5']
  data,ids = pickle.load(open(args.resultsFolderJoint+"results.p"))
  for i in range(len(zArr)):
    os.system("mkdir "+args.resultsFolderJoint+"sklearn_GMM_1_z_"+zStrArr[i])
    os.system("mkdir "+args.resultsFolderJoint+"sklearn_GMM_1_z_"+zStrArr[i]+"/sklearn_GMM_1_noLimits/")
    os.system("mkdir "+args.resultsFolderJoint+"sklearn_GMM_3_z_"+zStrArr[i])
    os.system("mkdir "+args.resultsFolderJoint+"sklearn_GMM_3_z_"+zStrArr[i]+"/sklearn_GMM_3_noLimits/")
    os.system("cp "+args.resultsFolderJoint+"results.p "+args.resultsFolderJoint+"sklearn_GMM_1_z_"+zStrArr[i]+"/")
    os.system("cp "+args.resultsFolderJoint+"results.p "+args.resultsFolderJoint+"sklearn_GMM_3_z_"+zStrArr[i]+"/")
    for j in range(len(ids)):
      if data['redshift']['value'][j] > zArr[i] and data['redshift']['value'][j] <= zArr[i]+deltaZ:
        if data['mass']['err'][0,j] < 0.8:
          os.system("cp "+args.resultsFolderJoint+"sklearn_GMM_1_noLimits/"+str(ids[j])+"_BEAGLE_GMM.p "+\
                    args.resultsFolderJoint+"sklearn_GMM_1_z_"+zStrArr[i]+"/sklearn_GMM_1_noLimits/")
          os.system("cp "+args.resultsFolderJoint+"sklearn_GMM_3_noLimits/"+str(ids[j])+"_BEAGLE_GMM.p "+\
                    args.resultsFolderJoint+"sklearn_GMM_3_z_"+zStrArr[i]+"/sklearn_GMM_3_noLimits/")


#Now perfrom Kelly07 fits
if args.runKelly07:
  if args.regularMlimGrid:
    massLimArr = np.array([7,7.5,8,8.5,9,9.5,10.0])
#    massLimArr = np.array([5])
   # massLimArr = np.array([9])
  elif args.regularMlimGridHigh:
    massLimArr = np.array([8.5,9,9.5,10,10.5,11])
  else:
    massLimits = pickle.load(open(args.resultsFolderJoint+"massLimits.p","r"))

    massLimArr = []
    for key in massLimits:
      if key not in massLimArr:
        massLimArr.append(massLimits[key]['massLimit'])

  inputArr = []
  for m in massLimArr:
    massStr = str(round(m,2))
    massStr = massStr.replace(".","p")
    plotSuffix = "_massLowerLim_"+massStr
    delayed = False
    if fileList:
      for i in range(len(folderListData['folder'])):
        #fileTest = os.path.isfile(folderListData['folder'][i]+"sklearn_GMM_1_noLimits/Gibbs_results"+plotSuffix+".p")
        #if fileTest == False:
        #  inputArr.append({'folder':folderListData['folder'][i],'massLim':str(m),'nGauss':'1'})
        fileTest = os.path.isfile(folderListData['folder'][i]+"sklearn_GMM_3_noLimits/Gibbs_results"+plotSuffix+".p")
        if fileTest == False:
          inputArr.append({'folder':folderListData['folder'][i],'massLim':str(m),'nGauss':'3', 'delayed':args.delayed, 'mTot':args.mTot})
    else:
      if args.regularZgrid:
        #for z in ['1p5','2p5','3p5','4p5','5p5']:
        for z in ['1p5']:
          fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_1_z_"+z+"/sklearn_GMM_1_noLimits/Gibbs_results"+plotSuffix+".p")
          #if fileTest == False:
          #  inputArr.append({'folder':args.resultsFolderJoint+"sklearn_GMM_1_z_"+z+"/",'massLim':str(m),'nGauss':'1','delayed':args.delayed})
          fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_3_z_"+z+"/sklearn_GMM_3_noLimits/Gibbs_results"+plotSuffix+".p")
          if fileTest == False:
            inputArr.append({'folder':args.resultsFolderJoint+"sklearn_GMM_3_z_"+z+"/",'massLim':str(m),'nGauss':'3','delayed':args.delayed, 'mTot':args.mTot})
      else:
        fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_1_noLimits/Gibbs_results"+plotSuffix+".p")
        if fileTest == False:
          inputArr.append({'folder':args.resultsFolderJoint,'massLim':str(m),'nGauss':'1', 'delayed':args.delayed, 'mTot':args.mTot})
        fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_3_noLimits/Gibbs_results"+plotSuffix+".p")
#        fileTest = False
#horrible Emma tweak here - should add re-run command line option
        if fileTest == False:
#        if 1==1:
          inputArr.append({'folder':args.resultsFolderJoint,'massLim':str(m),'nGauss':'3', 'delayed':args.delayed, 'mTot':args.mTot})
        #print 'before run_Kelly07', fileTest
        #run_Kelly07({'folder':args.resultsFolderJoint,'massLim':str(m),'nGauss':'3'})
        #print 'here'

  #print inputArr
  #run_Kelly07(inputArr[0])
  #sys.exit()
  #print 'inputArr: ', inputArr
  if len(inputArr) > 0:
    pool.map(run_Kelly07, inputArr)

#Also fit using David's algorithm
if args.runSimplifyFits:
  if fileList:
    inputArr = []
    for i in range(len(folderListData['folder'])):
      print "mkdir "+folderListData['folder'][i]+"fitsFiles/"
      os.system("mkdir "+folderListData['folder'][i]+"fitsFiles/")
      inputArr.append({'folder':folderListData['folder'][i],'sfrFree':args.sfrFree})

    pool.map(simplify_fits, inputArr)

  else:
    os.system("mkdir "+args.resultsFolderJoint+"fitsFiles/")
    command = "python -W ignore simplify_fits_files.py -r "+args.resultsFolderJoint+" -o "+args.resultsFolderJoint+"fitsFiles/"
    if not args.sfrFree:
      command = command + " --add-sfr"
    os.system(command)


if args.runDavids:
  if args.regularMlimGrid:
    massLimArr = np.array([7,7.5,8,8.5,9,9.5,10.0])
  else:
    massLimits = pickle.load(open(args.resultsFolderJoint+"massLimits.p","r"))
    massLimArr = []
    for key in massLimits:
      if key not in massLimArr:
        massLimArr.append(massLimits[key]['massLimit'])
  inputArr = []
  if fileList:
    for i in range(len(folderListData['folder'])):
      for m in massLimArr:
        inputArr.append({'folder':folderListData['folder'][i]+'fitsFiles/','massLimit':m,'sfrFree':args.sfrFree, 'mTot':args.mTot})
      os.system('cp '+folderListData['folder'][i]+"results.p "+folderListData['folder'][i]+"fitsFiles/")
    pool.map(run_modelling_massLims, inputArr)
  else:
    for m in massLimArr:
      inputArr.append({'folder':args.resultsFolderJoint+'fitsFiles/','massLimit':m,'sfrFree':args.sfrFree, 'mTot':args.mTot})
      #run_modelling_massLims({'folder':args.resultsFolderJoint+'fitsFiles/','massLimit':m,'sfrFree':args.sfrFree, 'mTot':args.mTot})
      #sys.exit()
    os.system('cp '+args.resultsFolderJoint+"results.p "+args.resultsFolderJoint+"fitsFiles/")
    print inputArr
    #run_modelling_massLims(inputArr[-1])
    pool.map(run_modelling_massLims, inputArr)


#and now provide measurements for Salmon, Santini and with medians
if args.runMeasureMedians:
  if args.regularMlimGrid:
    massLimArr = np.array([7,7.5,8,8.5,9,9.5,10.0])
  elif args.regularMlimGridHigh:
    massLimArr = np.array([8.5,9,9.5,10,10.5,11])
  else:
    massLimits = pickle.load(open(args.resultsFolderJoint+"massLimits.p","r"))
    massLimArr = []
    for key in massLimits:
      if key not in massLimArr:
        massLimArr.append(massLimits[key]['massLimit'])
  print massLimArr
  if fileList:
    inputArr = []
    inputPhotCatData = np.genfromtxt(args.inputPhotCat,dtype=None,names=True)
    for i in range(len(folderListData['folder'])):
      for m in massLimArr: # TODO: currently doesn't work properly for other mass limits
        inputArr.append({'folder':folderListData['folder'][i], 'inputPhotCat':inputPhotCatData['cat'][i],\
                         'filterInfo':args.filterInfo,\
                         'mTot':args.mTot, 'massLim':str(m), 'objectID':args.idFlip})
  else:
    inputArr = []
    for m in massLimArr: # TODO: currently doesn't work properly for other mass limits
        print 'm: ', m
        inputArr.append({'folder':args.resultsFolderJoint, 'inputPhotCat':args.inputPhotCat, \
                         'filterInfo':args.filterInfo,\
                         'mTot':args.mTot, 'massLim':str(m), 'objectID':args.idFlip})
        measure_medians(inputArr[-1])
  #measure_medians(inputArr[-1])
  #pool.map(measure_medians, inputArr)

#Now perfrom Kelly07 fits with sigma allowed to vary
if args.runKelly07sigma:
  if args.regularMlimGrid:
    massLimArr = np.array([7,7.5,8,8.5,9,9.5])
    #massLimArr = np.array([9.0])
  else:
    massLimits = pickle.load(open(args.resultsFolderJoint+"massLimits.p","r"))
    massLimArr = []
    for key in massLimits:
      if key not in massLimArr:
        massLimArr.append(massLimits[key]['massLimit'])
  inputArr = []
  for m in massLimArr:
    massStr = str(round(m,2))
    massStr = massStr.replace(".","p")
    plotSuffix = "_massLowerLim_"+massStr
    if fileList:
      for i in range(len(folderListData['folder'])):
        fileTest = os.path.isfile(folderListData['folder'][i]+"sklearn_GMM_1_noLimits/Gibbs_results_sigma"+plotSuffix+".p")
        if fileTest == False:
          inputArr.append({'folder':folderListData['folder'][i],'massLim':str(m),'nGauss':'1'})
        fileTest = os.path.isfile(folderListData['folder'][i]+"sklearn_GMM_3_noLimits/Gibbs_results_sigma"+plotSuffix+".p")
        if fileTest == False:
          inputArr.append({'folder':folderListData['folder'][i],'massLim':str(m),'nGauss':'3'})
    else:
      #fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_1_noLimits/Gibbs_results_sigma"+plotSuffix+".p")
      #if fileTest == False:
      #  inputArr.append({'folder':args.resultsFolderJoint,'massLim':str(m),'nGauss':'1'})
      fileTest = os.path.isfile(args.resultsFolderJoint+"sklearn_GMM_3_noLimits/Gibbs_results_sigma"+plotSuffix+".p")
      #if fileTest == False:
      if 'true' == 'true':
        inputArr.append({'folder':args.resultsFolderJoint,'massLim':str(m),'nGauss':'3'})


  pool.map(run_Kelly07_sigma_vary, inputArr)
