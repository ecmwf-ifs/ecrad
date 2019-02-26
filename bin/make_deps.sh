#!/bin/sh

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
