from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat

#setup
fermipath = "fermi/swiftj0243.fits.gz" # barytime v. amplitude for flux; barytime v. frequency for spin
swiftflux = "swift/SwiftJ0243.6p6124.lc.fits" # time
pathname = Path("/opt/data/mirror")

fermifilename = list(pathname.glob(fermipath))[0]
swiftfilename = list(pathname.glob(swiftflux))[0]
fdata = fits.open(fermifilename)
sdata, sheader = fits.getdata(swiftfilename,  header=True)
tzero = swiftbat.met2mjd(swiftbat.string2met('2017-09-04T00:00:00'))
print(tzero)
fdata = fdata[2].data
print(len(fdata['BARYTIME']))

fmidt = (fdata['BARYTIME'][:-1] + fdata['BARYTIME'][1:]) / 2
#i_f_s gives indices for fdata that match to sdata
i_f_s = np.searchsorted(sdata['TIME'], fmidt)
sdataflux = sdata['RATE'][i_f_s]
sdatatime = sdata['TIME'][i_f_s]

changeftime = np.diff(fdata['BARYTIME'])
changefreqs = np.diff(fdata['FREQUENCY']) / changeftime

# freq_diff 
freq_diff = np.searchsorted(fdata['BARYTIME'], fmidt)
fdataflux = fdata['AMPLITUDE'][freq_diff]
shiftfdataflux = fdataflux * 0.04133440433120353

print(len(fmidt))
print(len(sdataflux), len(fdataflux))
print(len(changefreqs), len(changeftime))

fig, axes = plt.subplots(nrows = 2, ncols = 1)
# plot Fermi Flux by Swift Flux
axes[0].plot(fdataflux, sdataflux, ".")
# find line of best fit
a, b = np.polyfit(fdataflux, sdataflux, 1)
axes[0].plot(fdataflux, a*fdataflux+b)  
axes[0].set(xlabel = 'Flux Fermi', ylabel = 'Flux Swift')
axes[0].legend([f"Slope: {a}"])

# plot freq vs Fermi flux and freq vs Swift flux
axes[1].plot(changefreqs, sdataflux, ".")
# axes[1].plot(changefreqs, fdataflux, ".")
axes[1].plot(changefreqs, shiftfdataflux, ".")
axes[1].set(xlabel = 'Rate of Change in Frequency', ylabel = 'Flux')
axes[1].legend(['Swift', 'Fermi * slope'])

plt.show()