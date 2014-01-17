#!/bin/bash
rm -f *.txt
export GSA_ROOT_DIR=/home/karthikg/broad/gsa-unstable
#-Djava.library.path is needed if you are using JNI_LOGLESS_CACHING, else not needed
java  -Djava.library.path=${GSA_ROOT_DIR}/PairHMM_JNI -jar ${GSA_ROOT_DIR}/dist/GenomeAnalysisTK.jar  -T HaplotypeCaller \
-R /data/broad/samples/joint_variant_calling/broad_reference/Homo_sapiens_assembly19.fasta \
-I /data/simulated/sim1M_pairs_final.bam \
--dbsnp /data/broad/samples/joint_variant_calling/dbSNP/00-All.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 10.0 \
--pair_hmm_implementation JNI_LOGLESS_CACHING \
-o output.raw.snps.indels.vcf

#--pair_hmm_implementation JNI_LOGLESS_CACHING \
#-I /data/simulated/sim1M_pairs_final.bam \
#-I /data/broad/samples/joint_variant_calling/NA12878_low_coverage_alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam \
