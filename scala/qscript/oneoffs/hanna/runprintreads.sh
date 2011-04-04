#!/bin/sh
java -Djava.io.tmpdir=/broad/shptmp/hanna -jar ~/src/Sting/dist/Queue.jar \
    --script PrintReadsAcrossManySamples.q \
    -gatk ~/src/Sting/dist/GenomeAnalysisTK.jar \
    -R /humgen/1kg/reference/human_g1k_v37.fasta \
    -I ~/tests/1600samples/1kg_t2d.list --max_bams 2000 --step_size 10 -bsub -jobQueue week $1
