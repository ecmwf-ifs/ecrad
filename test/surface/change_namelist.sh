#!/bin/bash
# Change entries in a Fortran namelist file
# Usage: change_namelist.sh infile.nam outfile.nam key1=value1 key2=value2 ...

INFILE=$1
OUTFILE=$2
shift
shift

SEDLINE=""
while [ "$1" ]
do
    FOUND=$(echo $1 | grep '=')
    if [ ! "$FOUND" ]
    then
	echo "Error in $0: argument '$1' not of the form key=value"
	exit 1
    fi
    KEY=$(echo $1 | awk -F= '{print $1}')
    VALUE=$(echo $1 | awk -F= '{print $2}')
    FOUND=$(grep $KEY $INFILE)
    if [ ! "$FOUND" ]
    then
	echo "Error: $KEY not found in $INFILE"
	exit 1
    fi

    SEDLINE="$SEDLINE -e s|^[[:space:]]*"$KEY".*|"$KEY"="$VALUE",|"
    shift
done
#echo sed $SEDLINE $INFILE ">" $OUTFILE
sed $SEDLINE $INFILE > $OUTFILE
