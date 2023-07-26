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

print(f'{pathname}*{mspattern}')
datafiles = sorted(glob.glob('/opt/data/mirror/swift/*/bat/rate/sw*brtms.lc.gz'))
# print(datafiles)
print(len(datafiles))
# print(list(pathname.joinpath(obsid).glob(mspattern)))

def prev_fast_FFT_len(n):
    ntry = n
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

# main list of data
alldata = []

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
    #for i in range(10,20):
        #print(frate[i],freqs[i])

    # component amplitude * 2 is generally how far away from the curve the data points are
    ncomponents = 20
    # sort from lowest to highest and then flip array
    ibest = np.argsort(np.abs(frate))[::-1]

    # Save FFT results into txt file
    for i in ibest[:ncomponents]:
        templist = [f"{i}", f"{freqs[i]}", f"{np.abs(frate[i])}", f"{np.rad2deg(np.angle(frate[i]))}"]
        # alldata.append(templist)
        file.write(f"{templist}")
        file.write("\n")
        # file.write(f"{i:7d} {freqs[i]:8.4f} {np.abs(frate[i]):10.4f} {np.rad2deg(np.angle(frate[i])):5.5f} \n")

file.close()

# Using for loop
lines = []
with open('FFT.txt', 'r') as f:
    for line in f:
        lines.append(line.strip())

# print(lines)

