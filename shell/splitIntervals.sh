#!/bin/sh

# Splits an interval list into multiple files

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input> <output> [ .. <output> ]"
    exit 1
else
    inputFile=$1
    shift

    totalLines=$(wc -l < $inputFile)

    exec 3< $inputFile

    numHeaders=0
    while true; do
	read -u 3 nextLine
	if [ $? -ne 0 ]; then break; fi
	if [[ $nextLine != @* ]]; then break; fi
	((numHeaders++))
    done

    numFiles=$#
    ((numIntervals = totalLines - numHeaders))

    if [ $numIntervals -lt $numFiles ]; then
	echo "Error: Number of intervals $numIntervals is less than the number of files $numFiles."
	exec 3<&-
	exit 1
    fi

    ((linesPerFile = numIntervals / numFiles))
    ((remainder = numIntervals % numFiles))

    ((linesPerFile++))

    fileNumber=0
    for outputFile in $@; do

	# Earlier files with get the remainder until it's no longer needed.
	if [ $fileNumber -eq $remainder ]; then ((linesPerFile--)); fi
	((fileNumber++))

	head -n $numHeaders $inputFile > $outputFile

	for ((line=0; line<$linesPerFile; line++)); do
	    echo "$nextLine" >> $outputFile
	    read -u 3 nextLine
	    if [ $? -ne 0 ]; then break; fi
	done
    done

    exec 3<&-
fi
