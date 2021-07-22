import numpy as np
import argparse
import os
import warnings
from pathos.multiprocessing import ProcessingPool
warnings.filterwarnings("ignore")

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
                    '--delayed',
                    help = "fit GMM and Kelly07 to instantaneous SFRs for a delayed history",
                    action="store_true",
                    default=False,
                    dest="delayed",
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
    command = 'python -W ignore run_GMM_2d_on_cluster.py -f '+args.resultsFolder+' --n-gauss 3'
    if args.mTot:
        command = command + ' --m-tot'
    if args.delayed:
        command = command + ' --delayed'
    print(command)
    os.system(command)

#perfrom Kelly07 fits
if args.runKelly07:
    command =   'python -W ignore run_GMM_2d_on_cluster.py -f '+args.resultsFolder+' --n-gauss 3 --n-chains 1 --min-iter 30 --max-iter 30 --n-K 3'
    if args.mTot:
        command = command + ' --m-tot'
    if args.delayed:
        command = command + ' --delayed'
    print(command)
    os.system(command)
