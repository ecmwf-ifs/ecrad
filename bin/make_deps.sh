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


EXT="parkind1.o|yomhook.o|yomcst.o|yomdyncore.o|yomlun.o|abor1.o|yomtag.o|mpl_module.o|yommp0.o"

while [ "$1" ]
do
    DEPS=$(egrep -i '^[ \t]*use' $1 | awk '-F[ ,]' '{print $2".o"}' | tr '[:upper:]' '[:lower:]' | egrep -v "$EXT" | tr '\n' ' ')
    if [ "$DEPS" ]
    then
	echo $1 | awk -F. '{print $1"'".o: $DEPS"'"}'
    fi
    shift
done
