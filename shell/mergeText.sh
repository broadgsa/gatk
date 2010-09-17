#!/bin/sh

# Merges a set of files, skipping over common headers.

if [ $# -lt 2 ]; then
    echo "Usage: $0 <output> <input> [ .. <input> ]"
    exit 1
elif [ $# -eq 2 ]; then
    cp $2 $1
else
    outputFile=$1
    shift

    test -e $outputFile && rm -f $outputFile

    exec 3< $1
    exec 4< $2

    startLine=1
    while true; do
	read -u 3 header1
	if [ $? -ne 0 ]; then break; fi
	read -u 4 header2
	if [ $? -ne 0 ]; then break; fi
	if [ "$header1" != "$header2" ]; then break; fi
	echo "$header1" >> $outputFile
	((startLine++))
    done
    
    exec 3<&-
    exec 4<&-
    
    for inputFile in $@; do
	tail -n +$startLine $inputFile >> $outputFile
    done
fi
