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
    ncomponents = 10
    # sort from lowest to highest and then flip array
    
    # Oddly different shapes??
    # print(frate)
    # print(freqs.shape)
    # print((datasegment['TIME']).shape)
    # Save FFT results into txt file

    # use x for unsorted, i for sorted
    for x, i in zip(range(ncomponents), ibest[:ncomponents]):
        # file.write(f"{i:7d} {freqs[i]:8.4f} {np.abs(frate[i]):10.4f} {np.rad2deg(np.angle(frate[i])):5.5f} \n")
        # print(freqs[i])
        file.write(f"{datasegment['TIME'][x] - tzero} {freqs[x]} {np.abs(frate[x])} {np.rad2deg(np.angle(frate[x]))} \n")

file.close()
timelist = []
freqlist = []


f = open('FFT.txt', 'r')  # We need to re-open the file
for line in f:
    line = line.strip()
    columns = line.split()
    time  = float(columns[0])
    frequency = float(columns[1])
    timelist.append(time)
    freqlist.append(frequency)

# print(timelist)
# print(freqlist)

fig, ax = plt.subplots(nrows = 1, ncols = 1)
ax.plot(timelist, freqlist, ".")
ax.set(xlabel = "TIME", ylabel = "FREQUENCIES")
fig.tight_layout()
plt.show()
