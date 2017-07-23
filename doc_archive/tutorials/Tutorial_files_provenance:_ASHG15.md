## Tutorial files provenance: ASHG15

http://gatkforums.broadinstitute.org/gatk/discussion/6760/tutorial-files-provenance-ashg15

<p>This document is intended to be a record of how the tutorial files were prepared for the AHSG 2015 hands-on workshop.</p>
<hr />
<h3>Reference genome</h3>
<p>This produces a 64 Mb file (uncompressed) which is small enough for our purposes, so we don't need to truncate it further, simplifying future data file preparations.</p>
<pre><code class="pre_md"># Extract just chromosome 20
samtools faidx /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta 20 &gt; human_g1k_b37_20.fasta

# Create the reference index
samtools faidx human_g1k_b37_20.fasta

# Create sequence dictionary
java -jar $PICARD CreateSequenceDictionary R=human_g1k_b37_20.fasta O=human_g1k_b37_20.dict

# Recap files
-rw-rw-r-- 1 vdauwera wga      164 Oct  1 14:56 human_g1k_b37_20.dict
-rw-rw-r-- 1 vdauwera wga 64075950 Oct  1 14:41 human_g1k_b37_20.fasta
-rw-rw-r-- 1 vdauwera wga       20 Oct  1 14:46 human_g1k_b37_20.fasta.fai</code class="pre_md"></pre>
<hr />
<h3>Sequence data</h3>
<p>We are using the 2nd generation CEU Trio of NA12878 and her husband and child in a WGS dataset produced at Broad with files names after the library preps, Solexa-xxxxxx.bam.</p>
<h4>1. Extract just chromosome 20:10M-20M bp and filter out chimeric pairs with -rf BadMate</h4>
<pre><code class="pre_md">java -jar $GATK -T PrintReads -R /path/to/bundle/current/b37/human_g1k_v37_decoy.fasta -I /path/to/Solexa-272221.bam -o NA12877_wgs_20_10M20M.bam -L 20:10000000-20000000 -rf BadMate 

java -jar $GATK -T PrintReads -R /path/to/bundle/current/b37/human_g1k_v37_decoy.fasta -I /path/to/Solexa-272222.bam -o NA12878_wgs_20_10M20M.bam -L 20:10000000-20000000 -rf BadMate 

java -jar $GATK -T PrintReads -R /path/to/bundle/current/b37/human_g1k_v37_decoy.fasta -I /path/to/Solexa-272228.bam -o NA12882_wgs_20_10M20M.bam -L 20:10000000-20000000 -rf BadMate 

# Recap files
-rw-rw-r-- 1 vdauwera wga     36240 Oct  2 11:55 NA12877_wgs_20_10M20M.bai
-rw-rw-r-- 1 vdauwera wga 512866085 Oct  2 11:55 NA12877_wgs_20_10M20M.bam
-rw-rw-r-- 1 vdauwera wga     36176 Oct  2 11:53 NA12878_wgs_20_10M20M.bai
-rw-rw-r-- 1 vdauwera wga 502282846 Oct  2 11:53 NA12878_wgs_20_10M20M.bam
-rw-rw-r-- 1 vdauwera wga     36464 Oct  2 12:00 NA12882_wgs_20_10M20M.bai
-rw-rw-r-- 1 vdauwera wga 505001668 Oct  2 12:00 NA12882_wgs_20_10M20M.bam</code class="pre_md"></pre>
<h4>2. Extract headers and edit manually to remove all contigs except 20 and sanitize internal filepaths</h4>
<pre><code class="pre_md">samtools view -H NA12877_wgs_20_10M20M.bam &gt; NA12877_header.txt

samtools view -H NA12878_wgs_20_10M20M.bam &gt; NA12878_header.txt

samtools view -H NA12882_wgs_20_10M20M.bam &gt; NA12882_header.txt</code class="pre_md"></pre>
<p>Manual editing is not represented here; basically just delete unwanted contig SQ lines and remove identifying info from internal filepaths.</p>
<h4>3. Flip BAM to SAM</h4>
<pre><code class="pre_md">java -jar $PICARD SamFormatConverter I=NA12877_wgs_20_10M20M.bam O=NA12877_wgs_20_10M20M.sam

java -jar $PICARD SamFormatConverter I=NA12878_wgs_20_10M20M.bam O=NA12878_wgs_20_10M20M.sam

java -jar $PICARD SamFormatConverter I=NA12882_wgs_20_10M20M.bam O=NA12882_wgs_20_10M20M.sam

#Recap files
-rw-rw-r-- 1 vdauwera wga 1694169101 Oct  2 12:28 NA12877_wgs_20_10M20M.sam
-rw-rw-r-- 1 vdauwera wga 1661483309 Oct  2 12:30 NA12878_wgs_20_10M20M.sam
-rw-rw-r-- 1 vdauwera wga 1696553456 Oct  2 12:31 NA12882_wgs_20_10M20M.sam</code class="pre_md"></pre>
<h4>4. Re-header the SAMs</h4>
<pre><code class="pre_md">java -jar $PICARD ReplaceSamHeader I=NA12877_wgs_20_10M20M.sam O=NA12877_wgs_20_10M20M_RH.sam HEADER=NA12877_header.txt

java -jar $PICARD ReplaceSamHeader I=NA12878_wgs_20_10M20M.sam O=NA12878_wgs_20_10M20M_RH.sam HEADER=NA12878_header.txt    

java -jar $PICARD ReplaceSamHeader I=NA12882_wgs_20_10M20M.sam O=NA12882_wgs_20_10M20M_RH.sam HEADER=NA12882_header.txt    

# Recap files
-rw-rw-r-- 1 vdauwera wga 1694153715 Oct  2 12:35 NA12877_wgs_20_10M20M_RH.sam
-rw-rw-r-- 1 vdauwera wga 1661467923 Oct  2 12:37 NA12878_wgs_20_10M20M_RH.sam
-rw-rw-r-- 1 vdauwera wga 1696538104 Oct  2 12:38 NA12882_wgs_20_10M20M_RH.sam</code class="pre_md"></pre>
<h4>5. Sanitize the SAMs to get rid of MATE_NOT_FOUND errors</h4>
<pre><code class="pre_md">java -jar $PICARD RevertSam I=NA12877_wgs_20_10M20M_RH.sam O=NA12877_wgs_20_10M20M_RS.sam SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=false REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false ATTRIBUTE_TO_CLEAR=null SANITIZE=true MAX_DISCARD_FRACTION=0.001

java -jar $PICARD RevertSam I=NA12878_wgs_20_10M20M_RH.sam O=NA12878_wgs_20_10M20M_RS.sam SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=false REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false ATTRIBUTE_TO_CLEAR=null SANITIZE=true MAX_DISCARD_FRACTION=0.001

java -jar $PICARD RevertSam I=NA12882_wgs_20_10M20M_RH.sam O=NA12882_wgs_20_10M20M_RS.sam SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=false REMOVE_DUPLICATE_INFORMATION=false REMOVE_ALIGNMENT_INFORMATION=false ATTRIBUTE_TO_CLEAR=null SANITIZE=true MAX_DISCARD_FRACTION=0.001

# Recap files
-rw-rw-r-- 1 vdauwera wga 1683827201 Oct  2 12:45 NA12877_wgs_20_10M20M_RS.sam
-rw-rw-r-- 1 vdauwera wga 1652093793 Oct  2 12:49 NA12878_wgs_20_10M20M_RS.sam
-rw-rw-r-- 1 vdauwera wga 1688143091 Oct  2 12:54 NA12882_wgs_20_10M20M_RS.sam</code class="pre_md"></pre>
<h4>6. Sort the SAMs, convert back to BAM and create index</h4>
<pre><code class="pre_md">java -jar $PICARD SortSam I=NA12877_wgs_20_10M20M_RS.sam O=NA12877_wgs_20_10M20M_V.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE

java -jar $PICARD SortSam I=NA12878_wgs_20_10M20M_RS.sam O=NA12878_wgs_20_10M20M_V.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE

java -jar $PICARD SortSam I=NA12882_wgs_20_10M20M_RS.sam O=NA12882_wgs_20_10M20M_V.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE

#recap files
-rw-rw-r-- 1 vdauwera wga     35616 Oct  2 13:08 NA12877_wgs_20_10M20M_V.bai
-rw-rw-r-- 1 vdauwera wga 508022682 Oct  2 13:08 NA12877_wgs_20_10M20M_V.bam
-rw-rw-r-- 1 vdauwera wga     35200 Oct  2 13:06 NA12878_wgs_20_10M20M_V.bai
-rw-rw-r-- 1 vdauwera wga 497742417 Oct  2 13:06 NA12878_wgs_20_10M20M_V.bam
-rw-rw-r-- 1 vdauwera wga     35632 Oct  2 13:04 NA12882_wgs_20_10M20M_V.bai
-rw-rw-r-- 1 vdauwera wga 500446729 Oct  2 13:04 NA12882_wgs_20_10M20M_V.bam</code class="pre_md"></pre>
<h4>7. Validate BAMs; should all output &quot;No errors found&quot;</h4>
<pre><code class="pre_md">java -jar $PICARD ValidateSamFile I=NA12877_wgs_20_10M20M_V.bam M=SUMMARY

java -jar $PICARD ValidateSamFile I=NA12878_wgs_20_10M20M_V.bam M=SUMMARY

java -jar $PICARD ValidateSamFile I=NA12882_wgs_20_10M20M_V.bam M=SUMMARY</code class="pre_md"></pre>