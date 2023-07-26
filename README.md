# ICR-research
Repository for the project from the Institude for Computing in Research --- Analysis of Swift J0243.6+6124

## Background

## Getting Started

### Packages
numpy, astropy, scipy, matplotlib, FTOOLS, swiftbat, swifttools, BatAnalysis

Most can be installed with `pip install`
BatAnalysis can be installed via https://github.com/parsotat/BatAnalysis
FTOOLS installation instructions can be found here: https://heasarc.gsfc.nasa.gov/docs/software/ftools/ftools_menu.html

Barycorr is also used for barycenter correction of the data, so you will need to set up the CALDB calibration system. Instructions can be found here:  https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/caldb_remote_access.html

### Data

To access the starting Swift data:
`wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2017_11//00059158012/`

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

To view plots for refining the frequency via rate v. phase relationship, do:
```
python3 freqsearch.py 
```