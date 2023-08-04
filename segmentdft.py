from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat

#setup
timecolname = 'BARYTIME'
mspattern ="j024_bary/*.bary" 
pathname = Path("/opt/data/mirror/swift")

# mspattern ="bat/rate/sw*brtms.lc.gz"
# pathname = Path("/opt/data/mirror/swift")
# obsid="00059158012"
# filename = list(pathname.joinpath(obsid).glob(mspattern))[0]

datafiles = sorted(pathname.glob(mspattern))
print("number of datafiles:", len(datafiles))
tb = 0.064    # size of the timebin (ignoring jumps)
tzero = swiftbat.string2met('2017-09-26T00:00:00')

cycles = 10
fapprox = 0.1018

Nsample = int(cycles / fapprox / tb)
print(Nsample)

# dft on data give dft value at frequency with phase reference to tzero
def dft(times, data, freq, tzero):
    theta = (((times-tzero) * freq) % 1) * 2 * np.pi # in units of radians from 0 to 2pi
    phasor = np.exp(-1j * theta)
    phasor -= phasor.mean()
    data = data-data.mean()
    # Discrete Fourier transform value (complex) at the freq
    signal = np.sum(phasor * data)/len(data)
    print(np.angle(signal)) # what does this do??
    return signal

midtimes = []
components = []
for f in datafiles:
    data, header = fits.getdata(f,  header=True)
    # splitting the data by time gaps
    splitlocs = np.argwhere(np.diff(data[timecolname]) > 1.5*tb).ravel() + 1 # gives index of the next time block

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

        # Do a DFT for each Nsample stretch of datasegment:

    
        df= -.12/600 -.9/12000 - .1/30000
        # -4.8e-3
        # -1.6e-3
        #-2.5e-4  # Adjustment goes here 

        for istart in np.arange(0, len(datasegment)-Nsample, Nsample):
            subsegment = datasegment[istart:istart+Nsample]
            rate = np.sum(subsegment['COUNTS'][:,0:2], axis=-1)/tb
            # tmid = subsegment[timecolname][len(subsegment)//2] - tzero
            # print(tmid, rate, tzero)
            component = dft(times = subsegment[timecolname], data=rate, tzero=tzero, freq=fapprox + df)
            midtimes.append(subsegment[timecolname][Nsample//2])
            components.append(component)
fig, axes = plt.subplots(nrows = 2, ncols = 1, sharex = True)
axes[0].plot(midtimes, np.unwrap(np.angle(components)) / (2 * np.pi), '.')

isort = np.argsort(midtimes)
midtimes = np.array(midtimes)[isort]
components = np.array(components)[isort]
        # Plot up the phase angle in cycles vs time
axes[1].plot(midtimes, np.unwrap(np.angle(components)) / (2 * np.pi), '.')

axes[0].set(xlabel = "Times (s)", ylabel = "Phase (Cycles)")
axes[0].set_title("Unsorted data")
axes[1].set(xlabel = "Times (s)", ylabel = "Phase (Cycles)")
# change in phase; horiz line = no change in freq + freq is correct, diag line = no change in freq + incorrect freq
# curved line = freq is changing + freq is correct where the tangent of the curved line is horiz (deriv = 0)
# increase slope toward the right = freq is higher
# slope (cycles/sec) = real freq - fapprox 
pass
# np.unwrap undoes jumps at 2pi=0 radians, so […6.1, 6.2, 0.02, 0.12…] becomes […6.1,6.2,6.3,6.4…]
plt.show()
# increasing phase = increasing spin frequency of star