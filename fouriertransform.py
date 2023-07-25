from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat


mspattern ="bat/rate/sw*brtms.lc.gz"
pathname = Path("/opt/data/mirror/swift")
obsid="00059158012"

filename = list(pathname.joinpath(obsid).glob(mspattern))[0]
data, header = fits.getdata(filename,  header=True)
tb = np.median(np.diff(data['TIME']))    # size of the timebin (ignoring jumps)
# print(header)
# print(data)
# print(tb)

# splitting the data by time gaps
splitlocs = np.argwhere(np.diff(data['TIME']) > 1.5*tb).ravel() + 1 # gives index of the next time block
allsegs = np.array([np.array(list()) for segments in np.split(data, splitlocs)])

for datasegment in np.split(data, splitlocs):
    starttime = swiftbat.met2datetime(datasegment['TIME'][0]) # time of spacecraft, not always accurate because of clock error
    duration = datasegment['TIME'].ptp()
    print(f"{starttime:%Y-%m-%dT%H:%M:%S} + {duration:5.0f} seconds to the end of the block")
    if duration > 1300:
    # if duration < 400:
        longdatasegment = datasegment
        # TODO: check if this is the longest data segment

allsegs = np.array(allsegs)        
        
# lower number for FFT
def prev_fast_FFT_len(n):
    ntry = n
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

n = prev_fast_FFT_len(len(longdatasegment) - int(60 / tb))
# Trim to the last n values (throwing away data that is probably during slew)
datasegment = longdatasegment[-n:]
duration = datasegment['TIME'].ptp()
print(f"{duration:.3f} seconds after trimming")

# Rate for the segment
rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb
'''
fig,ax = plt.subplots(nrows = 1, ncols = 1)
ax.plot(np.arange(len(rate)) * tb, rate)
ax.set(xlabel='Time??', ylabel='Rate (counts/sec)')
fig.tight_layout()
'''

# ABOUT FOURIER TRANSFORM
# Forward (time series->frequency) FFT for real input
# Subtract the mean to avoid a huge term at 0 frequency

# subtracting the mean makes it easier to see variations instead of just average
frate = sp.fft.rfft(rate - np.mean(rate), norm = "forward")
# print(frate[10:20])
# not transformed, used as reference for interpreting FFTs
freqs = sp.fft.rfftfreq(len(rate), tb)
# print(freqs[10:20])
for i in range(10,20):
    print(frate[i],freqs[i])

'''
fig,ax = plt.subplots(1,1)
ax.plot(freqs, np.abs(frate))
ax.set(xlabel="Frequency (Hz)", ylabel="Amplitude")
fig.tight_layout()
'''

# Forward (time series->frequency) FFT for real input
# Subtract the mean to avoid a huge term at 0 frequency
frate = sp.fft.rfft(rate, norm = "forward")
# print(frate[10:20])
freqs = sp.fft.rfftfreq(len(rate), tb)
# print(freqs[10:20])

invrate = sp.fft.irfft(frate, norm = "forward")
#for i in range(10,20):
#    print(frate[i],freqs[i])
times = np.arange(len(rate)) * tb


estrate = np.zeros(len(rate))
# component amplitude * 2 is generally how far away from the curve the data points are
ncomponents = 20
# sort from lowest to highest and then flip array
ibest = np.argsort(np.abs(frate))[::-1]
# count of cycles in entire observation

'''
fig,(axf,axr) = plt.subplots(nrows = 2, ncols = 1)
axf.plot(freqs, np.abs(frate))
axf.set(xlabel="Frequency (Hz)", ylabel="Amplitude")
axr.plot(times, rate, ".")
for i in ibest[:ncomponents]:
    print(f"{i:8d} {freqs[i]:8.4f} {np.abs(frate[i])} {np.rad2deg(np.angle(frate[i]))}")
    estrate += (1 if i == 0 else 2) * np.real(frate[i] * np.exp(times * freqs[i] * 1j * 2 * np.pi))
    # phase: angle of the spin of the star as of the start of the datasegment (only at correct frequency)
    axr.plot(times, estrate)
axr.set(xlabel="Seconds", ylabel="Rate", xlim = [0,20])
fig.tight_layout()
'''

tzero = swiftbat.string2met('2017-11-01T00:00:00')
print(tzero)
fapprox = 0.1019 # + divide .005 by seconds between each time chunk (1000 data points apart)
rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb
cycle, phase = np.divmod((datasegment['TIME']-tzero) * fapprox, 1) # spin angle of a specific point, seconds since tzero
fig, axes = plt.subplots(nrows = 2, ncols = 1, sharex = True)
print(len(phase), len(rate))

# plotting just the 1300s segment
for i in range(4):
    axes[0].plot(phase[1000*i:1000+1000*i], rate[1000*i:1000*i+1000], ".") # rate = brightness
    print(phase[1000*i])
axes[0].set(xlabel="Phase (in cycles)", ylabel="Rate (brightness)")
axes[0].set_title(f'{fapprox}')

# break the segment into 4 pieces
for segs in allsegs:
    segpieces = 4
    pointsperplot = 2 * int(1/(0.064 * fapprox))
    # for each segment, plot 2 cycles of data
    for istart in np.arange(0, len(segs), 1 + len(segs) // segpieces):
        sl = slice(istart, istart + pointsperplot)
        axes[1].plot(phase[sl], rate[sl], ".", label=f"{datasegment['time'][sl.start] - tzero:.0f}") 

axes[1].legend()
axes[1].set_title("")
    
fig.tight_layout()

plt.show()
