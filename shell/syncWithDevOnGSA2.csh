#!/bin/tcsh

setenv HERE "java tribble scala analysis"
setenv THERE \~/dev/GenomeAnalysisTKFromLaptop/trunk

rsync -e ssh -aCvz $HERE depristo@gsa1:$THERE
