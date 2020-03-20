#!/bin/bash

# This script sits in the bin directory of an ecRad distribution. To
# create an ecRad practical session in the present working directory,
# type the full path to this script.  It will then copy over the
# relevant files or create symbolic links.

# Stop if an error occurs
set -ex

ECRADDIR=$(dirname $0)/..
PRACDIR=$ECRADDIR/practical

ln -s -f $ECRADDIR/bin/ecrad
ln -s -f $ECRADDIR/data
ln -s -f $PRACDIR/era5slice.nc
ln -s -f $PRACDIR/ecrad_practical.pdf

cp -a -f $PRACDIR/ecradplot .
cp -f $PRACDIR/plot_*.py $PRACDIR/compare_*.py .
cp -f $PRACDIR/clean.sh $PRACDIR/config.nam $PRACDIR/README .

set +ex

echo 
echo "*** You are now ready to start the ecRad practical - congratulations! ***"
