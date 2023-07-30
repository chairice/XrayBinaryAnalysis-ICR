from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat
import glob

#setup
mspattern ="bat/rate/sw*brtms.lc.gz"
pathname = Path("/opt/data/mirror/swift")
file = open("FFT.txt","w")
# file.write("Data_Pt.  Freq    Amplitude Phase_Angle \n")
tzero = swiftbat.string2met('2017-09-01T00:00:00')

print(f'{pathname}*{mspattern}')
datafiles = sorted(glob.glob('/opt/data/mirror/swift/*/bat/rate/sw*brtms.lc.gz'))
# print(datafiles)
# print(len(datafiles))
# print(list(pathname.joinpath(obsid).glob(mspattern)))

def prev_fast_FFT_len(n):
    ntry = n
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

for f in datafiles:
    data, header = fits.getdata(f,  header=True)
    # size of the timebin (ignoring jumps)
    tb = np.median(np.diff(data['TIME']))
    
    # splitting the data by time gaps
    # gives index of the next time block
    splitlocs = np.argwhere(np.diff(data['TIME']) > 1.5*tb).ravel() + 1
    # print(f)
    # print out all the time gaps and the durations
    for datasegment in np.split(data, splitlocs):
        # time of spacecraft, not always accurate because of clock error
        starttime = swiftbat.met2datetime(datasegment['TIME'][0])
        n = prev_fast_FFT_len(len(datasegment) - int(60 / tb))
        datasegment = datasegment[-n:]
        duration = datasegment['TIME'].ptp()
        # print(f"{starttime:%Y-%m-%dT%H:%M:%S} + {duration:5.10f} s after trimming")

        # Rate for the segment
        rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb

        # subtracting the mean makes it easier to see variations instead of just average
        frate = sp.fft.rfft(rate - np.mean(rate), norm = "forward")
        # not transformed, used as reference for interpreting FFTs
        freqs = sp.fft.rfftfreq(len(rate), tb)
    
        # component amplitude * 2 is generally how far away from the curve the data points are
        ibest = np.argsort(np.abs(frate))[::-1]
        ncomponents = 1
        # sort from lowest to highest and then flip array
    
        # Save FFT results into txt file
        # the time at the midpoint because FFT returns n/2 frequencies
        tmid = datasegment['TIME'][len(datasegment)// 2] - tzero # how does moving this to the midddle change the accuracy of the data?
        for i in ibest[:ncomponents]:
            freq = freqs[i]
            if 0.2 < freq < 0.24:
                freq = freq/2
            if 0.1 < freq < 0.12:
                file.write(f"{tmid:7.3f} {freq:0.5f} {np.abs(frate[i]):3.4f} \n")
                break

file.close()
timelist = []
freqlist = []
amplist = []

# Reopen file to add time and freqs to a list
f = open('FFT.txt', 'r')
for line in f:
    line = line.strip()
    columns = line.split()
    time  = float(columns[0])
    frequency = float(columns[1])
    amps = float(columns[2])
    amplist.append(amps)
    timelist.append(time)
    freqlist.append(frequency)
    
'''
bestfreq = np.argsort(amplist)[::-1]

for x, i in zip(range(len(freqlist)), bestfreq):
    timelist[x] = timelist[bestfreq[i]]
    freqlist[x] = freqlist[bestfreq[i]]

print(timelist[:ncomponents])
print(freqlist[:ncomponents])
'''

fig, ax = plt.subplots(nrows = 1, ncols = 1)
ax.plot(timelist, freqlist, ".")
ax.set(xlabel = "TIME (seconds)", ylabel = "FREQUENCIES (Hz)")
fig.tight_layout()
plt.show()
