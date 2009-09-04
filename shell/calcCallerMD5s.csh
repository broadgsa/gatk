#!/bin/tcsh

foreach model (ONE_STATE THREE_STATE EMPIRICAL)
java -Xmx4096m -jar ~/dev/GenomeAnalysisTK/trunk/dist/GenomeAnalysisTK.jar -l INFO -L 1:10,000,000-11,000,000 -R /broad/1KG/reference/human_b36_both.fasta -T SingleSampleGenotyper -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -varout ./$model.geli.calls -lod 5 -m $model
end
wc -l *.geli.calls
md5sum *.geli.calls


