from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat

#setup
mspattern ="bat/rate/sw*brtms.lc.gz"
pathname = Path("/opt/data/mirror/swift")
obsid="00059158012"

filename = list(pathname.joinpath(obsid).glob(mspattern))[0]
data, header = fits.getdata(filename,  header=True)
tb = np.median(np.diff(data['TIME']))    # size of the timebin (ignoring jumps)

# splitting the data by time gaps
splitlocs = np.argwhere(np.diff(data['TIME']) > 1.5*tb).ravel() + 1 # gives index of the next time block

# print out all the time gaps and the durations
for datasegment in np.split(data, splitlocs):
    # time of spacecraft, not always accurate because of clock error
    starttime = swiftbat.met2datetime(datasegment['TIME'][0])
    duration = datasegment['TIME'].ptp()
    print(f"{starttime:%Y-%m-%dT%H:%M:%S} + {duration:5.10f} seconds to the end of the block")
    if duration > 1300:
        longdatasegment = datasegment      


# get the highest number of datapoints for FFT
def prev_fast_FFT_len(n):
    ntry = n
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

# refind the amount of datapoints afte cutting off the first minute (BAT is still adjusting the first min)
n = prev_fast_FFT_len(len(longdatasegment) - int(60 / tb))

# Trim to the last n values
datasegment = longdatasegment[-n:]
duration = datasegment['TIME'].ptp()
print(f"{duration:.3f} seconds after trimming")

# adjusting work
tzero = swiftbat.string2met('2017-11-01T00:00:00')
print(tzero)
fapprox0 = 0.1018
for df in np.arange(0, 1e-5, 0.25e-5):
    # fig, axes = plt.subplots(nrows = 11, ncols = 1, sharex = True)
    fig, axes = plt.subplots(nrows = 1, ncols = 1)
    # counts per second since tzero
    shiftpersecond = 1
    ax = axes
    fig.suptitle(f"{df:g}")
    rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb
    
    segpieces = 4
    # Break the segment into 4 pieces
    # 1/cycles/samples = samples/cycle
    ncyclesperplot = 4
    # ncyclesperplot * samples gives ncycles
    pointsperplot = ncyclesperplot * int(1/(0.064 * fapprox0))

    # 2 cycle of data for all 11 data segments
    for i, datasegment in enumerate(np.split(data, splitlocs)):
        datasegment = datasegment[len(datasegment)//2:]
        datasegment = datasegment[:pointsperplot]
        segrate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb
        segcycle, segphase = np.divmod((datasegment['TIME']-tzero) * (fapprox0 + df), 1)
        shift = (datasegment[0]['TIME'] - tzero) * shiftpersecond
        #y-axis units: rate but shifted by the start time since tzero
        # subtracting mean removes any potential background noise that may shift data downwards
        ax.plot(segphase, segrate + shift - np.mean(segrate), ".")
    
    fig.tight_layout()

plt.show()
