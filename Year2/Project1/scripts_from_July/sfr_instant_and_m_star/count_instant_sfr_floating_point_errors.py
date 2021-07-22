

import os
import numpy as np
from astropy.io import fits
from scipy import integrate

fields = ['0A2744C', '1A2744P', '2M0416C', '3M0416P', '4M0717C', '5M0717P', '6M1149C', '7M1149P']
fields = ['1A2744P']

for field in fields:
    print('')
    print(field)
    folder = '/home/ls861/rds/rds-rm665-beagle_shared/lester_shared/lester_results/{}/fit_001/'.format(field)
    
    folder = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_jul_2020/from_cluster/fields/1/'
    
    fileList = os.listdir(folder)
    for f, file in enumerate(fileList):
        if '.fits.gz' in file:
#            print(file, f, len(fileList))
            data = fits.open(folder+file)
            tau = np.power(10,data['POSTERIOR PDF'].data['tau'])
            age = np.power(10,data['POSTERIOR PDF'].data['max_stellar_age'])
   

            norm_denom2 = (-tau*np.exp(-age/tau)*(tau+age))  +np.power(tau,2)


            norm_denom = []
            
            for i in range(len(tau)):
                sfr_function = lambda t: t*np.exp(-t/tau[i])
    
                integral = integrate.quad(sfr_function, 0, age[i])
                norm_denom.append(integral[0])
#                print(norm_denom[i], norm_denom2[i])
                
                if integral[1]/integral[0] > 1e-5:
                    print('integration error', i, file, integral[0], integral[1]/integral[0])
                    
            norm_denom = np.array(norm_denom)


            
            print(' ')
            print(min(norm_denom2), max(norm_denom2))
            print(min(norm_denom), max(norm_denom))
            
            
            '''
            
            
            log10_age = np.log10(age)
            log10_mass = data['POSTERIOR PDF'].data['mass']            
            log10_norm = log10_mass - np.log10(norm_denom)

            exp_term = -age/tau
            exp_idx = (exp_term <= -100.0)
            exp_term[exp_idx] = -100.0
#            
            temp_sfr = log10_norm + log10_age + np.log10(np.exp(exp_term))
            sfr = np.where(temp_sfr < -30.0, -30.0, temp_sfr)

            '''

                    
                    
            data.close()


#            if len(temp_sfr[norm_idx]) > 0 or len(temp_sfr[exp_idx]) > 0:
#                print(field, file, len(temp_sfr), len(temp_sfr[norm_idx]), len(temp_sfr[exp_idx]))
#              

#print(np.log10(1e-250))         # -250 
#print(np.log10(np.exp(-500)))   # -217
                    

#print(np.log10(1e-100))         # -100                  
#print(np.log10(np.exp(-200)))   # -87                    
                    
#print(np.log10(1e-50))         # -50             
#print(np.log10(np.exp(-100)))   # -43   
            
            
        
            
            
            
            
            
