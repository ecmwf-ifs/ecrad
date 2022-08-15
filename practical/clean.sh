#!/bin/bash
#
# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


# This script cleans a directory that has been used for the ecRad
# practical

set -x

# Remove plots
rm -f era5slice*.png

# Remove edited config files
rm -f config_*.nam

# Remove ecRad output netCDF files, i.e. all netCDF files except for
# era5slice.nc
rm -f $(ls -1 *.nc | egrep -v "era5slice.nc|era5slice_hydromet.nc")

# Remove Python cached files
rm -rf ecradplot/__pycache__

# Remove emacs autosave files
rm -f *~

