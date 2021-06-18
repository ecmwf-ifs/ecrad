#!/bin/sh
#
# (C) Copyright 2015- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


while [ "$1" ]
do
    HEADERS=$(grep '^#include' $1 | awk '-F"' '{print $2}')
    for FILE in $HEADERS
    do
	touch ../include/$FILE
    done
    shift
done
