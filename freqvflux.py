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
fdata, fheader = fits.getdata(fermifilename,  header=True)
sdata, sheader = fits.getdata(swiftfilename,  header=True)
tzero = swiftbat.met2mjd(swiftbat.string2met('2017-09-04T00:00:00'))
print(tzero)

# print out all the time gaps and the durations
for datasegment in np.split(data, splitlocs):
    # time of spacecraft, not always accurate because of clock error
    starttime = swiftbat.met2datetime(datasegment[timecolname][0])
    duration = datasegment[timecolname].ptp()
    print(f"{starttime:%Y-%m-%dT%H:%M:%S} + {duration:5.10f} seconds to the end of the block")   

    # refind the amount of datapoints afte cutting off the first minute (BAT is still adjusting the first min)
    n = len(datasegment) - int(60 / tb)

    # Trim to the last n values
    datasegment = datasegment[-n:]
    duration = datasegment[timecolname].ptp()
    print(f"{duration:.3f} seconds after trimming")



    df= -.12/600 -.9/12000 - .1/30000
plt.show()