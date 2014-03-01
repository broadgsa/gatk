#!/bin/bash
rm -f *.txt *.log
GSA_ROOT_DIR=/home/karthikg/broad/gsa-unstable
pair_hmm_implementation="VECTOR_LOGLESS_CACHING";
if [ "$#" -ge 1 ];
then
  pair_hmm_implementation=$1;
fi

#-Djava.library.path is needed if you wish to override the default 'packed' library
#java -jar $GSA_ROOT_DIR/target/GenomeAnalysisTK.jar   -T HaplotypeCaller \
java  -Djava.library.path=${GSA_ROOT_DIR}/public/VectorPairHMM/src/main/c++ -jar $GSA_ROOT_DIR/target/GenomeAnalysisTK.jar   -T HaplotypeCaller \
--dbsnp /data/broad/samples/joint_variant_calling/dbSNP/00-All.vcf \
-R /opt/Genomics/ohsu/dnapipeline/humanrefgenome/human_g1k_v37.fasta \
-I /data/simulated/sim1M_pairs_final.bam \
-stand_call_conf 50.0 \
-stand_emit_conf 10.0 \
--pair_hmm_implementation $pair_hmm_implementation \
-o output.raw.snps.indels.vcf 

#--pair_hmm_implementation JNI_LOGLESS_CACHING \
#-XL unmapped	\
#-I /data/simulated/sim1M_pairs_final.bam \
#-I /data/broad/samples/joint_variant_calling/NA12878_low_coverage_alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam \
#-I /data/broad/samples/joint_variant_calling/NA12878_high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam \
#-R /data/broad/samples/joint_variant_calling/broad_reference/Homo_sapiens_assembly18.fasta \
#-R /data/broad/samples/joint_variant_calling/broad_reference/Homo_sapiens_assembly19.fasta \
#-R /data/broad/samples/joint_variant_calling/broad_reference/ucsc.hg19.fasta \
#-R /opt/Genomics/ohsu/dnapipeline/humanrefgenome/human_g1k_v37.fasta \
#-R /data/broad/samples/joint_variant_calling/broad_reference/human_g1k_v37_decoy.fasta \
#--dbsnp /data/broad/samples/joint_variant_calling/dbSNP/00-All.vcf \
#--dbsnp /data/broad/samples/joint_variant_calling/dbSNP/dbsnp_138.hg19.vcf \
