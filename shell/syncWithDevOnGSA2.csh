#!/bin/tcsh

setenv HERE java
setenv THERE \~/dev/GenomeAnalysisTKFromLaptop/trunk

rsync -e ssh -aCvz $HERE depristo@gsa2:$THERE
