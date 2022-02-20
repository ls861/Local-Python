import legac_utils as legac
import numpy as np
import pickle
import matplotlib.pyplot as plt
import traceback
import datetime
from os import getcwd, listdir, path

# with open('.pickle', 'rb') as input_file:
#     x = pickle.load(input_file)

# with open('.pickle', 'wb') as output_file:
#     pickle.dump(x, output_file)

#%%

# =============================================================================
# SAVE PICKLE FILES
# =============================================================================


# directory = '/Users/lester/Documents/GitHub/Local-Python/Year2/Project2/Code/myspec/'
directory = './myspec/'

filenames_complete = []
for filename in listdir(directory):
    if 'myspec' in filename and '.pickle' in filename:
        filenames_complete.append(filename[7:-7]+'.fits')

directory = './Files/'
filenames_todo = []
for filename in listdir(directory):
    if '_spec1d_' in filename and filename not in filenames_complete:
        filenames_todo.append(filename)
        
print(filenames_todo) 

del_attrs = ['templates', 'templates_rfft', 'A_eq_templ', 'b_eq_templ', 'A_ineq_templ', 'b_ineq_templ', 'A_ineq_kinem', 'b_ineq_kinem', 'A_eq_kinem', 'b_eq_kinem', 'matrix']


for i in range(len(filenames_todo)):
    
    myspec = legac.spec1d.from_filename(directory+filenames_todo[i])
    
    plt.plot(myspec.wave, myspec.spec)
    plt.show()
    
    try:
        
        myspec.ppxf_fit(plot=True, clean=False) # clean=True is better but takes a long time.

    except Exception as E: # code to run if there is an error
    
        with open('./Errors/{}.err'.format(filenames_todo[i][:-5]), 'a') as f:

            f.write(str(datetime.datetime.now()) + '\n')
            f.write(str(E.args) + '\n') # E.args is the message of the error
            f.write(traceback.format_exc() + '\n') # can also add traceback, need import traceback
    
    else: # only run if no exception was raised
    
        for attr in del_attrs:
        	delattr(myspec.pp, attr) 
        
        with open('./myspec/myspec_{}.pickle'.format(filenames_todo[i][:-5]), 'wb') as output_file:
            pickle.dump(myspec, output_file)
    
    finally: # this runs regardless of the above
        plt.show()

# print(myspec.pp.__dict__['matrix'])







