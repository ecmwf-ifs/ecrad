#!/bin/bash

# This script cleans a directory that has been used for the ecRad
# practical

set -x

# Remove plots
rm -f era5slice*.png

# Remove edited config files
rm -f config_*.nam

# Remove ecRad output netCDF files, i.e. all netCDF files except for
# era5slice.nc
rm -f $(ls -1 *.nc | egrep -v era5slice.nc)

# Remove Python cached files
rm -rf ecradplot/__pycache__

# Remove emacs autosave files
rm -f *~

