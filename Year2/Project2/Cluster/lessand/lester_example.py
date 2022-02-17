import legac_utils as legac
import matplotlib.pyplot as plt

myspec = legac.spec1d.from_filename('legac_M32_v3.11_spec1d_54311.fits')
myspec.ppxf_fit(plot=True, clean=False) # clean=True is better but takes a long time.
myspec.continuum_subtraction()

fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)

ax0.plot(myspec.wave, myspec.spec-myspec.gas_model, 'k-', alpha=0.5)
ax0.plot(myspec.wave, myspec.cont, 'r-', alpha=0.5)
ax1.plot(myspec.wave, myspec.contsub_spec, 'k-', alpha=0.5)
ax1.plot(myspec.wave, myspec.gas_model, 'r-', alpha=0.5)
