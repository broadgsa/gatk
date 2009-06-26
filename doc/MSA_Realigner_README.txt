Overview:
The Multiple Sequence Alignment (MSA) Realigner is designed to consume one or more BAM files
and to fix reads which are misaligned due to the presence of insertions/deletions (indels)
within or near them.  While this isn't the forum for a detailed explanation of why indel-
containing reads get misaligned, it is important to note that an artifact of these misalignments
is what look to be SNPs near the site of the indel (however, since they aren't really SNPs, I
often refer to them as columns of mismatches at a position in a pileup).  These particular false
positive SNPs usually occur in clusters (generally defined as 2 or more mismatch columns within
n base pairs, where n is usually less than or equal to 10).  It is often the case that an aligner
will detect the indel in some of the reads and will fail to detect it in others; that is because
the aligners don't use knowledge about the other reads mapping to the same location when placing an
individual read.  It is the realigner's job to use all of the reads mapping to a given location to
find a consensus indel which best explains the data and which minimizes entropy within the reads.

There are 3-4 major steps to the realignment process:
Step 1: Determining (small) suspicious intervals which are likely in need of realignment
Step 2: Merging the intervals
Step 3: Running the realigner over those intervals
Optional Step 4: Rebuild your original BAM with cleaned reads

A more detailed explanation follows.

-----

Step 1: Determining (small) suspicious intervals which are likely in need of realignment

There are several methods for finding these intervals, which can be used in conjunction with one
another or separately.

A. In the case that aligners do find some reads with indels in them, one would want to make sure
that the other indel-containing reads in the pileup are aligned correctly.  Note that when using
aligners which don't allow for gapped alignments (e.g. MAQ with single-end reads) this method is
not useful.

Usage:
java -jar dist/GenomeAnalysisTK.jar -I <input.bam> -R <ref.fasta> -T IndelIntervals
-L <regionsToCheck.txt> -S SILENT -o <intervalsOutput1.txt>

Optional Arguments:
--minIndelsPerInterval N [the minimum number of indels at a given position necessary for emission; default=1]

--allow454Reads [don't filter out 454 reads (which inherently have false indels); default=false]


B. Occasionally it is the case that you have a SNP call set for your file that you'd like to use
in searching for clustered SNPS (which are suspicious).  Note that the realigner works best with
an unfiltered SNP list if at all possible.  The following method outputs clustered SNP intervals. 

Usage:
java -jar dist/GenomeAnalysisTK.jar -R <ref.fasta> -T SNPClusters
-B dbsnp,dbsnp,<input.rod>,eval,1KGSNPs,<SNPlist.txt> -o <intervalsOutput2.txt>

Optional Arguments:
--windowSize N [mismatch columns are considered clustered when they occur no more than N bp apart; default=10]


C. When you do not have (or do not want to use) an available SNP call set, the following method
outputs intervals of clustered mismatching intervals.  Generally, one would use method B or
method C, but not both.

Usage:
java -jar dist/GenomeAnalysisTK.jar -I <input.bam> -R <ref.fasta> -T MismatchIntervals
-L <regionsToCheck.txt> -S SILENT -o <intervalsOutput3.txt>

Optional Arguments:
--windowSize N [mismatch columns are considered clustered when they occur no more than N bp apart; default=10]

--allow454Reads [don't filter out 454 reads (which inherently have false indels); default=false]

--mismatchFraction f [fraction of reads that need to mismatch for the position to be considered mismatching; default=0.15]
Note that this fraction should be adjusted based on your particular data set.  For DEEP coverage and/or
when looking for indels with low allele frequency, this number should be smaller.


Step 2: Merging the intervals

At this point, you need to combine any intervals files you have into a
master list; this is done by running the interval merger.

Usage: java -jar dist/GenomeAnalysisTK.jar -I <input.bam> -R <ref.fasta> -T IntervalMerger
--intervalsToMerge intervalsOutput1.txt [--intervalsToMerge intervalsOutput2.txt]
[--intervalsToMerge intervalsOutput3.txt] -o <mergedIntervalList.txt>

Optional Arguments:
--allow454Reads [don't filter out 454 reads (which inherently have false indels); default=false]

--maxIntervalSize [max size in bp of merged intervals that we'll pass to the realigner; default=500]


Step 3: Running the realigner over your intervals
Usage: java -jar dist/GenomeAnalysisTK.jar -I <input.bam> -R <ref.fasta> -T IntervalCleaner
-L mergedIntervalList.txt -S SILENT

Optional Arguments:
--allow454Reads [don't filter out 454 reads (which inherently have false indels); default=false]

--OutputCleaned <output.bam> [the output BAM file to emit the reads; by default it writes all reads - whether or 
not they were realigned - which at all overlap the input intervals (but not those outside the intervals)]
--OutputCleanedReadsOnly [when used with OutputCleaned it instructs the realigner to emit ONLY realigned reads]
--bam_compression N [when used with OutputCleaned it determines the BAM compression; default=5, recommended=1]

--OutputIndels <indels.txt> [the output file (text) for the indels found]

--LODThresholdForCleaning d [LOD threshold above which the realigner will proceed to realign; default=5.0]
This term is equivalent to "significance" - i.e. is the improvement significant enough to merit realignment?
Note that this number should be adjusted based on your particular data set.  For LOW coverage and/or
when looking for indels with low allele frequency, this number should be smaller.

--EntropyThreshold f [percentage of mismatching base quality scores at a position to be considered having high entropy; default=0.15]
This is similar to the argument in the MismatchIntervals method.  The point here is that the realigner
will only proceed with the realignment (even above the given threshold) if it minimizes entropy among
the reads (and doesn't simply push the mismatch column to another position).  This parameter is just
a heuristic and should be adjusted based on your particular data set.

--maxConsensuses N [max alternate consensuses to try (necessary to improve performance in deep coverage); default=30]
If you need to find the optimal solution regardless of running time, use a higher number.

--maxReadsForConsensuses N [max reads (chosen randomly) used for finding the potential alternate consensuses
(necessary to improve performance in deep coverage); default=120]
If you need to find the optimal solution regardless of running time, use a higher number.


Optional Step 4: Rebuild your original BAM with cleaned reads
If you want your cleaned read BAM to contain ALL of the original reads
too (regardless of whether they were cleaned or fell within one of the
target intervals), you can do so in this last optional step.
Important note: this option works best with the
"-OutputCleanedReadsOnly" option in Step 3.

First, be sure to index your cleaned output BAM from Step 3:
samtools index <output.bam>

Usage: java -jar dist/GenomeAnalysisTK.jar -I <originalInput.bam> -R <ref.fasta> -T CleanedReadInjector
--cleaned_intervals mergedIntervalList.txt --cleaned_reads <previousOutput.bam> -S SILENT --output_bam <fullOutput.bam>

Optional Arguments:
--bam_compression N [when used with OutputCleaned it determines the BAM compression; default=5, recommended=1]


Questions or comments:
Email Eric Banks - ebanks@broadinstitute.org


