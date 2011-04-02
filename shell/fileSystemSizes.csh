#!/bin/tcsh

setenv DATE `date +"%m_%d_%Y"`
setenv RESULTS "fs_sizes.$DATE.txt"

rm -f $RESULTS
foreach fs ( /humgen/gsa-scr1/ /humgen/gsa-hphome1/ /humgen/gsa-hpprojects /humgen/gsa-lpprojects )
  du -sh $fs/* >> $RESULTS
end
