
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.integrate import quad

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def sfr_calc(sfh, mtot, msa, tau, tau_exp, alpha, beta, z):
    
    if sfh == 'DE':
        
        m       = 10**mtot                                                          # solar masses
        t       = 10**msa                                                           # yrs
        ta      = 10**tau                                                           # yrs

        A = m / ( -ta*np.exp(-t/ta)*(t+ta) + ta**2 )
        # https://www.wolframalpha.com/input/?i=A*%28t%29*exp%28-%28t%29%2FC%29
      
        sfr = A * t * np.exp(-t/ta)

    elif sfh == 'iDE':
        
        m       = 10**mtot                                                          # solar masses
        t       = 10**msa                                                           # yrs
        ta      = 10**tau                                                           # yrs
        A       = np.empty(len(m))
        
        for i in range(len(A)):            
            integrand = lambda T: T * np.exp(-T/ta[i])
            A[i] = m[i] / quad(integrand, 0, t[i])[0]
                
        sfr = A * t * np.exp(-t/ta)
           
    elif sfh == 'LE':
        
        m       = 10**mtot                                                          # solar masses
        t       = 10**msa                                                           # yrs
        ta      = 10**tau                                                           # yrs
        ta_exp  = 10**tau_exp                                                       # yrs
        A       = np.empty(len(m))
        
        for i in range(len(A)):          
            integrand = lambda T: (T*np.heaviside(ta[i] - T, 0)  + ta[i]*np.exp((ta[i]-T)/ta_exp[i])*np.heaviside(T - ta[i], 1))
            A[i] = m[i] / quad(integrand, 0, t[i])[0]
                
        sfr = A * (t*np.heaviside(ta - t, 0)  + ta*np.exp((ta-t)/ta_exp)*np.heaviside(t - ta, 1))
        
    elif sfh == 'DPL':
        
        m       = 10**mtot                                                          # solar masses
        ta      = 10**tau                                                           # yrs
        t = cosmo.age(z).value * 1E9                                                # yrs
        A       = np.empty(len(m))
        err     = np.empty(len(m))
        
        for i in range(len(A)):          
            integrand = lambda T: ( ( ((T/ta[i])**alpha[i]) + ((T/ta[i])**-beta[i]) ) ** -1 )
            A[i] = m[i] / quad(integrand, 0, t[i])[0]
            err[i] = quad(integrand, 0, t[i])[1] / quad(integrand, 0, t[i])[0]
        sfr = A * ( ( ((t/ta)**alpha) + ((t/ta)**-beta) ) ** -1 )
        
    else:
        sfr='ERROR in sfr_calc'
     
    return sfr, err


size = 15
fsize = 6

# =============================================================================
# get "real" parameters
# =============================================================================

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/fit/mock_catalogue_005_010.fits'
data_fits = fits.open(fileName)


mtot_r = np.log10(data_fits['GALAXY PROPERTIES'].data['m_tot'])                    # mtot is the value we selected as prior
#r_mstar = data_fits['GALAXY PROPERTIES'].data['m_star']
sfr_r = data_fits['STAR FORMATION'].data['SFR']


data_fits.close()


plt.figure(figsize=(fsize, fsize))
plt.title('Plot showing SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.scatter(mtot_r, np.log10(sfr_r))
plt.show()





# =============================================================================
# get BEAGLE parameters
# =============================================================================

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_006/astrodeep_002/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

#fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_007/astrodeep_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

fileName = '/Users/lester/Documents/GitHub/Local-Python/Astrodeep_feb_2020/from_cluster/param_008/astrodeep_001/pyp-beagle/data/BEAGLE_summary_catalogue.fits'

data_fits = fits.open(fileName)

id_b = np.asarray(data_fits[1].data['ID'], dtype=int) - 1
z_b = data_fits[1].data['redshift_mean']
z_med_b = data_fits[1].data['redshift_median']
mtot_b = data_fits[1].data['mass_mean']
#msa_b = data_fits[1].data['max_stellar_age_mean']
tau_b = data_fits[1].data['tau_mean']
alpha_b = data_fits[1].data['dpl_alpha_mean']
beta_b = data_fits[1].data['dpl_beta_mean']

data_fits.close()


# =============================================================================
# plot Redshifts
# =============================================================================

plt.figure(figsize=(fsize, fsize/2))
plt.title('DPL - Histogram of BEAGLE fitted redshifts', size=size)
plt.xlabel(r'Redshift', size=size)
plt.ylabel(r'Count', size=size)
plt.hist(z_b, bins=50, histtype=u'step', label='mean')
plt.hist(z_med_b, bins=50, histtype=u'step', label='median')
plt.legend()
plt.show()



# =============================================================================
# plot main sequence
# =============================================================================

sfr_b, err = sfr_calc('DPL', mtot_b, 0, tau_b, 0, alpha_b, beta_b, z_b)

plt.hist(err, bins=50)
plt.ylim(0, 10)

plt.show()



plt.figure(figsize=(fsize, fsize))
plt.title('DPL - Plot showing BEAGLE fitted SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
plt.scatter(mtot_r[id_b], np.log10(sfr_r[id_b]), label='input')
plt.scatter(mtot_b, np.log10(sfr_b), label='BEAGLE output')
plt.scatter(np.delete(mtot_r,id_b), np.log10(np.delete(sfr_r,id_b)), label='not fitted')
plt.legend()
plt.show()


plt.figure(figsize=(fsize, fsize))
plt.title('DPL - Plot showing BEAGLE fitted SFR vs Mass', size=size)
plt.xlabel(r'$\text{log}(m_{tot}/M_{\odot})$', size=size)
plt.ylabel(r'$\text{log}(\Psi / M_{\odot} yr^{-1})$', size=size)
plt.xlim(7.5, 11)
plt.ylim(-1, 3.5)
for i in range(len(sfr_b)):
    color = next(plt.gca()._get_lines.prop_cycler)['color']
    plt.plot( (mtot_r[id_b][i],mtot_b[i]), (np.log10(sfr_r[id_b][i]),np.log10(sfr_b[i])), '-', color=color)
    plt.plot( mtot_r[id_b][i], np.log10(sfr_r[id_b][i]), 'o', color=color, markersize=3)
plt.show()   
    

# This creates a list of indices for which the relative integration error is less than the specified value
ind = err <= 0.01
ind = [i for i, x in enumerate(ind) if x]









'''





# =============================================================================
# Resizing Arrays
# =============================================================================

z_AD = z_AD[id_B-1]
H160 = H160[id_B-1]
errH160 = errH160[id_B-1]
SN = abs(H160/errH160)

ind = (SN > 0) & (SN < 20)

z_B = z_B[ind]
z_Berr_m = z_Berr_m[ind]
z_Berr_p = z_Berr_p[ind]
z_AD = z_AD[ind]
H160 = H160[ind]
errH160 = errH160[ind]
SN = SN[ind]


# =============================================================================
# PLOT
# =============================================================================

pt = 20

#fig, ax = plt.subplots(figsize=np.array([6.4, 4.8]))
fig, ax = plt.subplots(figsize=np.array([17.6, 4.8]))
#fig.suptitle('BEAGLE vs Astrodeep redshifts', size=20, ha='center')
ax.set_title('BEAGLE vs Astrodeep redshifts', size=pt, ha='center')

h = ax.scatter(z_AD, z_B, c=SN, cmap='plasma', zorder=1)
ax.errorbar(z_AD, z_B, yerr=[z_Berr_m, z_Berr_p], linestyle="None", elinewidth=1, color='k', zorder=0)

x_lin = np.linspace(0, 20, 2)
ax.plot(x_lin, x_lin)

ax.set_xlim(0, 11)
ax.set_ylim(0, 13)

ax.set_xlabel(r'Astrodeep Redshift', size=pt)
ax.set_ylabel(r'BEAGLE Redshift', size=pt)

cbar = plt.colorbar(h)
cbar.set_label(label='SNR', size=pt)
cbar.ax.tick_params(labelsize=pt)
plt.show()

fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(SN, bins=50)


'''