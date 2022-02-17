import legac_utils as legac
import matplotlib.pyplot as plt
from os import getcwd, listdir

# myspec = legac.spec1d.from_filename('legac_M32_v3.11_spec1d_54311.fits')
# myspec.ppxf_fit(plot=True, clean=False) # clean=True is better but takes a long time.
# myspec.continuum_subtraction()

# fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)

# ax0.plot(myspec.wave, myspec.spec-myspec.gas_model, 'k-', alpha=0.5)
# ax0.plot(myspec.wave, myspec.cont, 'r-', alpha=0.5)
# ax1.plot(myspec.wave, myspec.contsub_spec, 'k-', alpha=0.5)
# ax1.plot(myspec.wave, myspec.gas_model, 'r-', alpha=0.5)


# print(myspec.pp.__dict__.keys())
# print(myspec.pp.gas_names)
# print(myspec.pp.gas_flux)
# print(myspec.pp.gas_flux_error)
# print(myspec.pp.velscale)


#%%
# for key in myspec.pp.__dict__.keys():
#     if key in ['fraction', 'ftol', 'gas_reddening', 'sky', 'reddening', 'reg_dim', 'tied', 'reddening_func', 'polyval', 'polyvander', 'A_eq_templ', 'b_eq_templ', 'A_ineq_kinem', 'b_ineq_kinem', 'A_eq_kinem', 'b_eq_kinem', 'gas_mpoly', 'chi2', 'polyweights', 'apoly'] or type(myspec.pp.__dict__[key]) in [int, bool, np.int64]:
#         pass
#     elif len(myspec.pp.__dict__[key])==6:
#         print(key, myspec.pp.__dict__[key])


#%%

import legac_utils as legac
import numpy as np
from os import getcwd, listdir, path
import pickle


if path.exists('errorfiles.pickle'):
    with open('errorfiles.pickle', 'rb') as input_file:
        errorfiles = pickle.load(input_file)
else:
    errorfiles = ['legac_M101_v3.11_spec1d_124654.fits', 
                  'legac_M101_v3.11_spec1d_129553.fits',
                  'legac_M101_v3.11_spec1d_130074.fits',
                  'legac_M101_v3.11_spec1d_131353.fits']

# legac_M101_v3.11_spec1d_125127.fits is extra lines error


if path.exists('gas_filenames.pickle'):
    with open('gas_filenames.pickle', 'rb') as input_file:
        gas_filenames = pickle.load(input_file)
    with open('gas_flux.pickle', 'rb') as input_file:
        gas_flux = pickle.load(input_file)
    with open('gas_flux_error.pickle', 'rb') as input_file:
        gas_flux_error = pickle.load(input_file)

    gas_names = list(gas_flux.keys())

    with open('gas_filenames_old.pickle', 'wb') as output_file:
        pickle.dump(gas_filenames, output_file)
    
    with open('gas_flux_old.pickle', 'wb') as output_file:
        pickle.dump(gas_flux, output_file)
        
    with open('gas_flux_error_old.pickle', 'wb') as output_file:
        pickle.dump(gas_flux_error, output_file)
        
else:
    gas_names = ['Hdelta', 'Hgamma', 'Hbeta', '[OII]3726', '[OII]3729', '[OIII]5007_d']
    gas_filenames = []
    gas_flux = {}
    gas_flux_error = {}
    
    for name in gas_names:
        gas_flux[name] = []
        gas_flux_error[name] = []
   
objects = 10

# for filename in listdir(getcwd()):

new_filenames = []
directory = '/export/data/fdeugenio/legac_team_share/spectra/'
# directory = './'
for filename in listdir(directory):
    if 'spec1d' in filename and filename not in gas_filenames and filename not in errorfiles:
        new_filenames.append(filename)

for filename in new_filenames[:objects]:
    print('LESTER STARTING: ', filename)
    errorfiles.append(filename)
    with open('errorfiles.pickle', 'wb') as output_file:
        pickle.dump(errorfiles, output_file)
        
    gas_filenames.append(filename)
    myspec = legac.spec1d.from_filename(directory+filename)
    myspec.ppxf_fit(plot=True, clean=False) # clean=True is better but takes a long time.
    
    gas_names_old = gas_names
    for name in [item for item in myspec.pp.__dict__['gas_names'] if item not in gas_names_old]:
        gas_names.append(name)
        gas_flux[name] = [-1.0 for i in range(len(gas_filenames)-1)]
        gas_flux_error[name] = [-1.0 for i in range(len(gas_filenames)-1)]    
    
    for i, name in enumerate(myspec.pp.__dict__['gas_names']):
        gas_flux[name].append(myspec.pp.__dict__['gas_flux'][i])
        gas_flux_error[name].append(myspec.pp.__dict__['gas_flux_error'][i])
    for name in list(set(gas_names).difference(myspec.pp.__dict__['gas_names'])):
        gas_flux[name].append(-1.0)
        gas_flux_error[name].append(-1.0)

    with open('gas_filenames.pickle', 'wb') as output_file:
        pickle.dump(gas_filenames, output_file)
    
    with open('gas_flux.pickle', 'wb') as output_file:
        pickle.dump(gas_flux, output_file)
        
    with open('gas_flux_error.pickle', 'wb') as output_file:
        pickle.dump(gas_flux_error, output_file)
        
    print('LESTER COMPLETE: ', filename)

    errorfiles = errorfiles[:-1]
    with open('errorfiles.pickle', 'wb') as output_file:
        pickle.dump(errorfiles, output_file)

'''
import pickle
with open('gas_filenames.pickle', 'rb') as input_file:
    gas_filenames = pickle.load(input_file)
    
print(gas_filenames)
print(len(gas_filenames))
'''
# with open('gas_flux.pickle', 'rb') as input_file:
#     gas_flux = pickle.load(input_file)
# print(gas_flux)



with open('gas_flux.pickle', 'rb') as input_file:
    gas_flux = pickle.load(input_file)


gas_names = list(gas_flux.keys())
print(gas_names)
print(gas_flux['Hbeta'])
print(gas_flux['[SII]6716'])
    
with open('errorfiles.pickle', 'rb') as input_file:
    errorfiles = pickle.load(input_file)
    
    
    
    
    
    
    
    
    
    
    