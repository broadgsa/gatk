#!/bin/bash


chr=$1

# set basic input and output paths. User should modify the following to point to:
# a) Input vcf,
# b) Output base path for all files
# c) Location of beagle.jar
# d) Path to reference file
beagle_prefix="CEUTSI.chr${chr}.bgl"
beagledir="/broad/shptmp/username/beagle/CEUTSI_Pilot1/chr${chr}/"
beaglejar="../../beagle/beagle.jar"

input_vcf="/broad/shptmp/username/beagle/CEUTSI_Pilot1/CEUTSI.recal.filtered.vcf"
tmpdir="/broad/shptmp/username/tmp"
bgloutprefix="recal.filtered"
reffile="/broad/1KG/reference/human_g1k_v37.fasta"


# Set to one to determine which sections of beagle pipeline to run
runinput=1 # Run VCF to Beagle input converter
runbgl=1 # Run Beagle
runoutput=1 # run Beagle output-to-VCF converter
runvarianteval=1 # Run Variant Evaluator to measure Beagle performance.

# Reference files for variant evaluator
dohapmap=1
do1kgchip=1

# Path to HapMap/1KG Chip truth sets for variant evaluator
if [ $dohapmap == 1 ]
then
    cmpfileh="/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/HapMap/3.2/genotypes_r27_nr.hg19_fwd.vcf"
fi

if [ $do1kgchip == 1 ]
then
    cmpfile1="/humgen/gsa-hpprojects/GATK/data/Comparisons/Validated/1kg_chip_jan2010/1000Genome.chip.hg19.filtered.vcf"
fi

outputbglvcf=$beagledir$bgloutprefix.$beagle_prefix.output.vcf


# check if Beagle directory exists. If not, create it.
if [ ! -d $beagledir ] 
then
    echo "Creating Beagle directory $beagledir"
    mkdir $beagledir
fi

if [ $runinput == 1 ] 
then
    echo "Running GATK to create Beagle input"

    java -Xmx4000m -jar ../dist/GenomeAnalysisTK.jar -L $chr -l INFO \
	-R $reffile -T ProduceBeagleInput \
	-B vcf,VCF,$input_vcf \
	-beagle $beagledir$beagle_prefix
fi

# now, run beagle
if [ $runbgl == 1 ]
then
    echo "Running Beagle..."
    java -Xmx8000m -Djava.io.tmpdir=$tmpdir -jar $beaglejar like=$beagledir$beagle_prefix out=$bgloutprefix

    # move output files to beagle directory
    cp ./$bgloutprefix.* $beagledir
    # unzip gzip'd files, force overwrite if existing
    gunzip -f $beagledir$bgloutprefix.$beagle_prefix.gprobs.gz
    gunzip -f $beagledir$bgloutprefix.$beagle_prefix.phased.gz
    #rename also Beagle likelihood file to mantain consistency
    cp $beagledir$beagle_prefix $beagledir$bgloutprefix.$beagle_prefix.like 
    cp $beagledir$beagle_prefix.log $beagledir$bgloutprefix.$beagle_prefix.log
fi

# run GATK to parse Beagle files and to produce output vcf
if [ $runoutput == 1 ]
then
    java -Xmx4000m -Djava.io.tmpdir=$tmpdir -jar ../dist/GenomeAnalysisTK.jar \
	-R $reffile -T BeagleOutputToVCF  -l INFO -L $chr \
	-B inputvcf,VCF,$input_vcf \
	-B beagleR2,BEAGLE,$beagledir$bgloutprefix.$beagle_prefix.r2 \
	-B beagleProbs,BEAGLE,$beagledir$bgloutprefix.$beagle_prefix.gprobs \
	-B beaglePhased,BEAGLE,$beagledir$bgloutprefix.$beagle_prefix.phased \
	-output $beagledir$bgloutprefix.$beagle_prefix.output.vcf 
fi

if [ $runvarianteval == 1 ]
then
    # finally, run VariantEval to produce useful comparisons between pre-and post-Beagle vcf's.
    if [ $dohapmap == 1 ]
    then
    java -Xmx4096m -jar ../dist/GenomeAnalysisTK.jar -T VariantEval \
	-R $reffile -l INFO -L $chr \
	-B eval_prebeagle,VCF,$input_vcf \
	-B eval_beagle,VCF,$outputbglvcf \
	-B comp_hapmap,VCF,$cmpfileh \
	-reportType Grep -o ${beagledir}$bgloutprefix.$beagle_prefix.variant_eval_hapmap_grep.txt
    fi

    if [ $do1kgchip == 1 ]
    then
    java -Xmx4096m -jar ../dist/GenomeAnalysisTK.jar -T VariantEval \
	-R $reffile -l INFO -L $chr \
	-B eval_prebeagle,VCF,$input_vcf \
	-B eval_beagle,VCF,$outputbglvcf \
	-B comp_1kgchip,VCF,$cmpfile1 \
	-reportType Grep -o ${beagledir}$bgloutprefix.$beagle_prefix.variant_eval_1kgchip_grep.txt
    fi
fi


