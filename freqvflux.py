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

print(len(sdata['RATE']))
print(len(fdata['AMPLITUDE']))

#i_f_s gives indices for fdata that match to sdata
i_f_s = np.searchsorted(sdata['TIME'], fdata['BARYTIME'])
sdata['RATE'] = sdata['RATE'][i_f_s]
sdata['TIME'] = sdata['TIME'][i_f_s]
# fdata['BARYTIME'][i_f_s]
# sdata['TIME']
                        
ratechangefreqs = []
ftimes = []

for i in np.arange(0, len(fdata['BARYTIME']) - 1):
    # get rate of change in frequency
    changetime = fdata['BARYTIME'][i+1] - fdata['BARYTIME'][i]
    changefreq = fdata['FREQUENCY'][i+1] - fdata['FREQUENCY'][i]
    ratefreqchange = changefreq / changetime
    ratechangefreqs.append(ratefreqchange)
    ftimes.append(changetime)

# ffluxoverfreq = []
# sfluxoverfreq = []

# for i in np.arange(0,249):
#     ffluxoverfreq.append(ratechangefflux[i] / ratechangefreqs[i]) 
#     sfluxoverfreq.append(ratechangesflux[i] / ratechangefreqs[i])


fig, axes = plt.subplots(nrows = 2, ncols = 1)
axes[0].plot(fdata['AMPLITUDE'], sdata['RATE'], ".")
#find line of best fit
a, b = np.polyfit(fdata['AMPLITUDE'], sdata['RATE'], 1)
axes[0].plot(fdata['AMPLITUDE'], a*fdata['AMPLITUDE']+b)  
print(a)
axes[0].set(xlabel = 'Flux Fermi', ylabel = 'Flux Swift')
axes[0].legend(['Slope: a'])

axes[1].plot(ratechangefreqs, ratechangefflux, ".")
axes[1].plot(ratechangefreqs, ratechangesflux, ".")
axes[1].set(ylabel = 'Change in Flux over Change in Frequency', xlabel = 'Time')

plt.show()