#!/bin/tcsh

setenv CMD "java -Djava.io.tmpdir=/broad/shptmp/depristo/tmp -jar /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/Queue.jar -statusTo depristo -S /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/analysis/depristo/distributedGATK/distributedGATKPerformance.scala -bsub --gatkjarfile /humgen/gsa-scr1/depristo/dev/GenomeAnalysisTKFromLaptop/trunk/dist/GenomeAnalysisTK.jar -dataset HiSeq $argv[2-$#argv]"

if ( $1 == 1 ) then
  pushd short; $CMD -jobQueue hour -run &
else if ( $1 == 2 ) then
  pushd long; $CMD -jobQueue gsa -long -run &
else
  $CMD
endif
