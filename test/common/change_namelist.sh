#!/bin/bash
#
# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Change entries in a Fortran namelist file
# Usage: change_namelist.sh infile.nam outfile.nam key1=value1 key2=value2 ...

INFILE=$1
OUTFILE=$2
shift
shift

SEDLINE=""
while [ "$1" ]
do
    if [[ "$1" != *"="* ]]
    then
        echo "Error in $0: argument '$1' not of the form key=value"
        exit 1
    fi
    KEY="${1%%=*}"
    VALUE="${1#*=}"
    if ! grep -q "$KEY" "$INFILE"
    then
        echo "Error: $KEY not found in $INFILE"
        exit 1
    fi

    SEDLINE="$SEDLINE -e s/^[[:space:]]*${KEY}.*/${KEY}=${VALUE},/"
    shift
done

echo "sed $SEDLINE \"$INFILE\" > \"$OUTFILE\""
sed $SEDLINE "$INFILE" > "$OUTFILE"
