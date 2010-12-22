#!/bin/tcsh

# what's the status of all of the gsa hosts
echo "GSA host status"
bhosts gsahosts

echo "\nGSA queue usage"
bjobs -u all -q gsa | awk '$2 !~ "USER" {print $2}' | sort | uniq -c

echo "\nGeneral computing resources"
bqueues gsa week short broad

echo "\nFH jobs"
bjobs -u gsa-adm

echo "\nFile system status"
ls /humgen/gsa-scr1 /humgen/1kg /humgen/gsa-hpprojects /humgen/gsa-hphome1 /humgen/gsa-pipeline /humgen/gsa-firehose2 /humgen/gsa-lpprojects

df -h /humgen/* /broad/shptmp 
