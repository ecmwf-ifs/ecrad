#!/bin/bash -e
#
# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Duplicate profiles with solar zenith angles varying from 0 to 90 degrees
# You need to have the nco netcdf tools in your PATH

INPUT=$1
OUTPUT=$2

NSZA=46
COS_SZA='1.0,0.999391,0.997564,0.994522,0.990268,0.984808,0.978148,0.970296,0.961262,0.951057,0.939693,0.927184,0.913545,0.898794,0.882948,0.866025,0.848048,0.829038,0.809017,0.788011,0.766044,0.743145,0.71934,0.694658,0.669131,0.642788,0.615661,0.587785,0.559193,0.529919,0.5,0.469472,0.438371,0.406737,0.374607,0.34202,0.309017,0.275637,0.241922,0.207912,0.173648,0.139173,0.104528,0.0697565,0.0348995,0.01'
#NSZA=2
#COS_SZA='1,0.1'

# Check for existence of NCO commands
command -v ncks >/dev/null 2>&1 || { \
 echo "###########################################################" ; \
 echo "### Error: NCO commands (ncks etc) needed but not found ###" ; \
 echo "###########################################################" ; \
 exit 1; }

ncks -O --mk_rec_dmn column $INPUT tmp0.nc
ncrcat -O -n $NSZA,1,0 tmp0.nc tmp_$OUTPUT
ncap2 -O -s "cos_solar_zenith_angle(:)={$COS_SZA}" tmp_$OUTPUT $OUTPUT
rm tmp0.nc tmp_$OUTPUT
