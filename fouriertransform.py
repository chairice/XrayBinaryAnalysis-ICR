from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat


mspattern ="bat/rate/sw*brtms.lc.gz"
pathname = Path("/opt/data/mirror/swift")
obsid="00059158012"

filename = list(pathname.joinpath(obsid).glob(mspattern))[0] # why zero??
data, header = fits.getdata(filename,  header=True)
tb = np.median(np.diff(data['TIME']))    # size of the timebin (ignoring jumps)
data, header = fits.getdata(filename, header = True)
# print(header)
# print(data)
print(dt)

# splitting the data by time gaps
splitlocs = 
