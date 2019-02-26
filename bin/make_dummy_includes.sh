#!/bin/sh

while [ "$1" ]
do
    HEADERS=$(grep '^#include' $1 | awk '-F"' '{print $2}')
    for FILE in $HEADERS
    do
	touch ../include/$FILE
    done
    shift
done
