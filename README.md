# ICR-research
Repository for the project from the Institude for Computing in Research --- Analysis of Swift J0243.6+6124

## Background
...

## Getting Started

### Packages
It's recommended to install Anaconda and set up a separate environment for this project. This is where you can install most of the required tools for this project. To do this, copy and paste the following into a new terminal:
- `conda create -n {environment name} python`
- `conda activate {environment name}`
- `pip install numpy, astropy, scipy, matplotlib, swiftbat, swifttools`

BatAnalysis can be installed via https://github.com/parsotat/BatAnalysis. 
FTOOLS installation instructions can be found here: https://heasarc.gsfc.nasa.gov/docs/software/ftools/ftools_menu.html 

Barycorr is also used for barycenter correction of the data, so you will need to set up the CALDB calibration system. Instructions can be found here:  https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_remote_access.html

### Data

To access the starting Swift data: 
```
wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2017_11//00059158012/
```

For additional data when refining freq, do:
```
python3 swiftJ0246_lightcurve.py
```

### Running the code

To view the plots for the data files in terms of rate or counts as a function of time, do:
```
python3 plottingTycho_200_200.py
```

To view plots for the Fourier transform of data and rate v. phase relationship, do:
```
python3 fouriertransform.py 
```

To find data on refined frequency via Discrete Fourier Transforms of barycentered data, do: 
```
python3 daybaryfreq.py
```
then 
```
cat barydayfreq.txt
```
or use another text editor to view the txt file. 