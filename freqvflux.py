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
pass

ratechangefreqs = []
ratechangesflux = []
ratechangefflux = []
ftimes = []
fmidtimes = []
# smidtimes = []

for i in np.arange(0, len(fdata['BARYTIME']) - 1):
    # get rate of change in frequency
    changeftime = fdata['BARYTIME'][i+1] - fdata['BARYTIME'][i]
    changeffreq = fdata['FREQUENCY'][i+1] - fdata['FREQUENCY'][i]
    ratefchange = changeffreq / changeftime
    ratechangefreqs.append(ratefchange)

    changestime = sdata['TIME'][i+1] - sdata['TIME'][i]
    changesfreq = sdata['RATE'][i+1] - sdata['RATE'][i]
    rateschange = changesfreq / changestime
    ratechangesflux.append(rateschange)

    # get rate of change in flux
    changeflux = fdata['AMPLITUDE'][i+1] - fdata['AMPLITUDE'][i]
    rateflux = changeflux / changeftime
    ratechangefflux.append(rateflux)
    ftimes.append(changeftime)
    fmidtimes.append(fdata['BARYTIME'][i] + changeftime / 2)
    # smidtimes.append(sdata['TIME'][i] + changestime / 2)

print(len(smidtimes))
pass

ffluxoverfreq = []
sfluxoverfreq = []

for i in np.arange(0,249):
    ffluxoverfreq.append(ratechangefflux[i] / ratechangefreqs[i]) 
    sfluxoverfreq.append(ratechangesflux[i] / ratechangefreqs[i])


fig, axes = plt.subplots(nrows = 2, ncols = 1)
axes[0].plot(ratechangefreqs, ratechangefflux, ".")
axes[0].plot(ratechangefreqs, ratechangesflux, ".")
axes[1].plot(fmidtimes, ffluxoverfreq, ".")
axes[1].plot(fmidtimes, sfluxoverfreq, ".")

axes[0].set(xlabel = 'Change in Frequency', ylabel = 'Change in Flux')
axes[0].legend(['fdata', 'sdata'])
axes[1].set(ylabel = 'Change in Flux over Change in Frequency', xlabel = 'Time')

plt.show()