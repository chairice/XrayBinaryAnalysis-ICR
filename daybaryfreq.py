from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat
import glob

#setup
timecolname = 'BARYTIME'
mspattern ="j024_bary/*.bary" 
pathname = Path("/opt/data/mirror/swift")
freqfile = open("/home/chairice/IXrayBinaryAnalysis-ICR/barydayfreq.txt","wt")

tzero = swiftbat.string2met('2017-09-26T00:00:00')
tfirst = tzero
tlast = swiftbat.string2met('2017-11-11T14:29:02')
tstarts = []

def prev_fast_FFT_len(n):
    ntry = n
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

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

#TODO: check if baryfreq and dayfreq are different, do dft on specific frequencies??
datafiles = sorted(pathname.glob(mspattern))
print("number of datafiles:", len(datafiles))
tb = 0.064
npointsday = 1024*1024 

for tzeroday in np.arange(tfirst, tlast, npointsday * tb / 2):
    dayrates = np.zeros(npointsday)

    for f in datafiles:
        data, header = fits.getdata(f,  header=True)    
        # gives index of the next time block
        # print(data.dtype)
        splitlocs = np.argwhere(np.diff(data[timecolname]) > 1.5*tb).ravel() + 1
        # print out all the time gaps and the durations
        for datasegment in np.split(data, splitlocs):
            # time of spacecraft, not always accurate because of clock error
            starttime = swiftbat.met2datetime(datasegment[timecolname][0])
            tstarts.append(starttime)
            timeinday = datasegment[timecolname][0] - tzeroday
            istart = int(timeinday / tb)
            if istart < 0 or istart > npointsday:
                continue

            if istart + len(datasegment) >= npointsday:
                continue

            # Rate for the segment
            rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/tb
            rate -= rate.mean()
            dayrates[istart:istart+len(datasegment)] = rate
            # print("day standard dev =", dayrates.stddev())

    # print(f"{min(tstarts)}, {max(tstarts)}")
    rate = dayrates

    # subtracting the mean makes it easier to see variations instead of just average
    frate = sp.fft.rfft(rate - np.mean(rate), norm = "forward")
    # not transformed, used as reference for interpreting FFTs
    freqs = sp.fft.rfftfreq(len(rate), tb)
            
    # component amplitude * 2 is generally how far away from the curve the data points are
    ibest = np.argsort(np.abs(frate))[::-1]
    ncomponents = 20
    # sort from lowest to highest and then flip array

    # Save FFT results into txt file
    # the time at the midpoint because FFT returns n/2 frequencies
    tmid = tzeroday + npointsday / 2 * tb
    for i in ibest:
        freq = freqs[i]
        # print(f"{tmid:7.3f} {freq:0.5f} {np.abs(frate[i]):3.4f} \n")
        if 0.2 < freq < 0.24:
            freq = freq/2
        if 0.1 < freq < 0.12:
            freqmessage = f"{tmid:7.3f} {freq:0.5f} {np.abs(frate[i]):3.4f}"
            print(freqmessage)
            freqfile.write(freqmessage+"\n")
            # print(freqmessage, file=freqfile)
            freqfile.flush()
            break



freqfile.close()