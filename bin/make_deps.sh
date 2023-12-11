#!/bin/bash
#
# (C) Copyright 2015- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set -e

EXT="___"

set -- $(echo $@ | tr -s ' ' '\n' | sort | xargs)
FILES=$(echo $@ | tr -s ' ' '\n' | sed 's/\.F90$/.o/g' | sort)
while [ "$1" ]
do
    DEPS=$(egrep -i '^[ \t]*use' $1 | sed 's/^ *//g' | awk '-F[ ,]' '{print $2".o"}' | tr '[:upper:]' '[:lower:]' | egrep -v "$EXT" | sort)
    #echo "Checking $1"
    #echo "files: $FILES" | xargs
    #echo "deps: $DEPS" | xargs
    DEPS=$(join <(echo "$DEPS") <(echo "$FILES") | sort -n | uniq | xargs)
    #echo "filtered: $DEPS" | xargs
    if [ "$DEPS" ]; then
        echo $1 | awk -F. '{print $1"'".o: $DEPS"'"}'
        #echo ""
        #echo ""
    fi
    shift
done
