import numpy as np
import argparse
import os
import warnings
from pathos.multiprocessing import ProcessingPool
warnings.filterwarnings("ignore")

def run_Kelly07(inputDict):
    print('HERE')
    success = False
    counter = 1
    while success == False and counter < 5:
        counter = counter+1
        command =   "python -W ignore run_GMM_2d_and_Kelly_sampler_cluster.py -f "+inputDict['folder']+\
                    " -r --n-gauss "+inputDict['nGauss']+" --no-limits --mass-lower "+inputDict['massLim'] +\
                    " --n-chains 1 --min-iter 30 --max-iter 30 --n-K 3"

        if inputDict['mTot']:
            command = command + ' --m-tot'

        print(command)
        try:
            os.system(command)
            success = True
        except:
            print('run_Kelly07 failed', counter)

parser = argparse.ArgumentParser()
parser.add_argument(
                    '-f', '--results-folder',
                    help="name of folder containing BEAGLE fits, give list of files in .txt file if more than one to be used to derive the mass limits (e.g. with the test scenarios)",
                    action="store",
                    type=str,
                    dest="resultsFolder",
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
                    '--regular-mlim-grid',
                    help = "set to true if you want a regular mlim grid",
                    action="store_true",
                    default=False,
                    dest="regularMlimGrid",
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

args = parser.parse_args()
pool = ProcessingPool(nodes=1)

#perform GMM fits
if args.runGMMfits:
    command = "python -W ignore run_GMM_2d_and_Kelly_sampler_cluster.py -f "+args.resultsFolder+" -g --n-gauss 3 --no-limits"
    if args.mTot:
        command = command + " --m-tot"
    print(command)
    os.system(command)

#perfrom Kelly07 fits
if args.runKelly07:
    if args.regularMlimGrid:
        massLimArr = np.array([1,7,7.5,8,8.5,9,9.5,10.0])

    inputArr = []
    for m in massLimArr:
        inputArr.append({'folder':args.resultsFolder,'massLim':str(m),'nGauss':'3', 'mTot':args.mTot})

    if len(inputArr) > 0:
        print(inputArr)
        pool.map(run_Kelly07, inputArr)
