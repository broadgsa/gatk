#!/usr/bin/csh
# ##############################################################################
# this file can be used to diff the flagstat values for a whole genome,
# along with sub-diffs for each chromosome.  This may be helpful if you're 
# debugging a specific reads problem, and you want to know the differences
# between the output of the GATK and samtools
# ##############################################################################

if ($#argv != 3) then
    echo "usage readStatsDiff <bamfile.bam> <referenceSequence.fasta> <outputFile>"; exit 1
endif

if (! -f $1) then
    echo bam file $2 does not seem to exist; exit 1
endif
if (! -f $2) then
    echo reference file $2 does not seem to exist; exit 1
endif
if ( -e $3) then
    echo output file $3 already exists!; exit 1
endif


# set the bam and the reference index
set bam = $1
set seq = $2
set out = $3

# diff the whole chromosome
java -ea -Xmx4096m \
-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=8015 \
-jar dist/GenomeAnalysisTK.jar \
-S SILENT \
-T FlagStatRead \
-I $bam \
-R $seq \
-oe 1.out 
samtools flagstat $bam > 2.out
echo "\n" >> $out
echo results for whole bam file >> $out
echo "\n" >> $out
diff 1.out 2.out >> $out
rm 1.out
rm 2.out



foreach chromosome (chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
java -ea -Xmx4096m \
-Xdebug -Xrunjdwp:transport=dt_socket,server=y,suspend=n,address=8015 \
-jar dist/GenomeAnalysisTK.jar \
-S SILENT \
-T FlagStatRead \
-I $bam \
-R $seq \
-oe 1.out \
-L $chromosome
samtools view -h -b $bam $chromosome > temp.bam
samtools flagstat temp.bam > 2.out
echo results for $chromosome >> $out
echo "\n" >> $out
diff 1.out 2.out >> $out
rm temp.bam
rm 1.out
rm 2.out
end
