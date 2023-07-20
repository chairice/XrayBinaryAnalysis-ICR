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
data, header = fits.getdata(filename, header = True)
# print(header)
# print(data)
print(tb)

# splitting the data by time gaps
splitlocs = np.argwhere(np.diff(data['TIME']) > 1.5*tb).ravel() + 1 # gives index of the next time block

for split in np.split(data, splitlocs):
    print(f"{swiftbat.met2datetime(split['TIME'][0]):%Y-%m-%dT%H:%M:%S} + {split['TIME'].ptp():5.0f} seconds to the end of the block")
