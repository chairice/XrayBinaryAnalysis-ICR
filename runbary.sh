#! /bin/env bash

thisdir=$(dirname $0)

topdir="/opt/data/mirror/swift"
barydir=${topdir}/j024_bary
mkdir $barydir

ra=40.9184386956300
dec=+61.4343770713400 

# orbitfile=${thisdir}/sw05000401003sao.fits

for origfile in ${topdir}/*/bat/rate/*brtms.lc.gz
do
    outfile=${barydir}/$(basename $origfile).bary
    orbitfile=$(dirname ${origfile})/../../auxil/*sao.fits.gz
    orbitfile=$(echo ${orbitfile})
    # echo ${orbitfile}
    # echo test${orbitfile}
    barycorr ra=${ra} dec=${dec} clockfile=CALDB infile=${origfile} outfile=${outfile} orbitfile=${orbitfile} barytime=yes
    # exit
done