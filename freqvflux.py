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

'''
for i in range(len(fdataflux)):
    if fdataflux[i] > 7:
        print(i)
'''

pass


fig, axes = plt.subplots(nrows = 2, ncols = 1)
# plot Fermi Flux by Swift Flux
axes[0].plot(fdataflux, sdataflux, ".")
# find line of best fit
a, b = np.polyfit(fdataflux[:19], sdataflux[:19], 1)
c, d = np.polyfit(fdataflux[19:64], sdataflux[19:64], 1)
axes[0].plot(fdataflux, a*fdataflux+b)  
axes[0].plot(fdataflux, c*fdataflux+d)  
axes[0].set(xlabel = 'Flux Fermi', ylabel = 'Flux Swift')
axes[0].legend([f"Slope: {a},{c}"])

# fdataflux = fdataflux * 0.04133440433120353

for flux in fdataflux:
    if flux < 7:
        flux = flux * a
    else:
        flux = flux * c
A = []
B = []
C = []
D = []
E = []

# get indices of <58065, 58140, 58460, 58600, then 60000
for i in np.arange(0, len(fmidt)):
    if fmidt[i] < 58065:
        A.append(i)
    if fmidt[i] >= 58065 and fmidt[i] < 58140:
        B.append(i)
    if fmidt[i] >= 58140 and fmidt[i] < 58460:
        C.append(i)
    if fmidt[i] >= 58460 and fmidt[i] < 58600:
        D.append(i)
    if fmidt[i] > 60000:
        E.append(i)

# plot freq vs Fermi flux and freq vs Swift flux
axes[1].plot(changefreqs, sdataflux, ".")
axes[1].plot(changefreqs[A], fdataflux[A], ".")
axes[1].plot(changefreqs[B], fdataflux[B], ".")
axes[1].plot(changefreqs[C], fdataflux[C], ".")
axes[1].plot(changefreqs[D], fdataflux[D], ".")
axes[1].plot(changefreqs[E], fdataflux[E], ".")
axes[1].set(xlabel = 'Rate of Change in Frequency', ylabel = 'Flux')
axes[1].legend(['Swift', 'Fermi * slope'])

# plot Fermi data
fig, axes = plt.subplots(nrows = 2, ncols = 1)
axes[0].plot(fmidt[A], changefreqs[A], ".")
axes[0].plot(fmidt[B], changefreqs[B], ".")
axes[0].plot(fmidt[C], changefreqs[C], ".")
axes[0].plot(fmidt[D], changefreqs[D], ".")
axes[0].plot(fmidt[E], changefreqs[E], ".")

axes[1].plot(fmidt[A], fdataflux[A], ".")
axes[1].plot(fmidt[B], fdataflux[B], ".")
axes[1].plot(fmidt[C], fdataflux[C], ".")
axes[1].plot(fmidt[D], fdataflux[D], ".")
axes[1].plot(fmidt[E], fdataflux[E], ".")

plt.show()