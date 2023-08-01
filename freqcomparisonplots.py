import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

listfiles = ['dayfreq.txt','barydayfreq.txt']
fig, axes = plt.subplots(nrows = 1, ncols = 1)
ax = axes

# plot FFT for all data with barycentered data and original data

for file in listfiles:
    timelist = []
    freqlist = []
    amplist = []

    # Reopen file to add time and freqs to a list
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        columns = line.split()
        time  = float(columns[0])
        frequency = float(columns[1])
        amps = float(columns[2])
        amplist.append(amps)
        timelist.append(time)
        freqlist.append(frequency)
        

    bestfreq = np.argsort(amplist)[::-1]

    for x, i in zip(range(len(freqlist)), bestfreq):
        timelist[x] = timelist[bestfreq[i]]
        freqlist[x] = freqlist[bestfreq[i]]

    ax.plot(timelist, freqlist, ".")
    ax.set(xlabel = "TIME (seconds)", ylabel = "FREQUENCIES (Hz)")
    
fig.tight_layout()
plt.legend(listfiles)
plt.show()

# axes[1].plot(freqlist, amplist, ".")
# axes[1].set(xlabel = "FREQUENCIES (Hz)", ylabel = "AMPLITUDE (|DFT of rate|)")