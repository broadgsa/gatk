## [How to] Generate a BAM for variant discovery (long)

http://gatkforums.broadinstitute.org/gatk/discussion/5969/how-to-generate-a-bam-for-variant-discovery-long

<h3>This document is an archived rough draft of <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6483">Tutorial#6483</a>. Please use the public tutorial. If you are interested in aligning to GRCh38, then please refer to a separate tutorial, <a href="https://software.broadinstitute.org/gatk/documentation/article?id=8017">Tutorial#8017</a>.</h3>
<hr />
<p>[work in progress--I am breaking this up into smaller chunks]
<a name="top"></a>
This document in part replaces the previous post <a href="http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format">(howto) Revert a BAM file to FastQ format</a> that uses HTSlib commands. The workflow assumes familiarity with the concepts given in <a href="http://gatkforums.broadinstitute.org/discussion/1317/collected-faqs-about-bam-files">Collected FAQs about BAM files</a>.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/72/6abbf031529f1a8287a302adb454aa.png" height="270"align="right" border="9"/>
<p>We outline steps to preprocess Illumina and similar tech DNA sequence reads for use in GATK's variant discovery workflow. This preprocessing workflow involves marking adapter sequences using <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters">MarkIlluminaAdapters</a> so they contribute minimally to alignments, alignment using the <a href="http://bio-bwa.sourceforge.net/bwa.shtml#3">BWA</a> aligner's maximal exact match (MEM) algorithm, and preserving and adjusting read and read meta data using <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a> for consistency and comparability of downstream results with analyses from the Broad Institute. With the exception of BWA, we use the most current versions of tools as of this writing. The workflow results in an aligned BAM file with appropriate meta information that is ready for processing with MarkDuplicates.</p>
<p>This workflow applies to three common types of sequence read files: (A) aligned BAMs that need realignment, (B) FASTQ format data and (C) raw sequencing data in BAM format. If you have raw data in BAM format (C), given appropriate read group fields, you can start with step 2. The other two formats require conversion to unmapped BAM (uBAM). We use Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#RevertSam">RevertSam</a> to convert an aligned BAM (A) or Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam">FastqToSam</a> to convert a FASTQ (B) to the uBAM. </p>
<p>We address options relevant to process reads extracted from an interval as well as options to process large files, in our case a ~150G file called <code>Solexa-272222</code>. The tutorial uses a smaller file of reads aligning to a genomic interval, called <code>snippet</code> derived from <code>Solexa-272222</code>, for faster processing. The example commands apply to the larger file. Some comments on the workflow: </p>
<ul>
<li>The workflow reflects a <em>lossless</em> operating procedure that retains original FASTQ read information within the final BAM file such that data is amenable to reversion and analysis by different means. These practices make scaling up and longterm storage efficient, as one needs only store the final BAM file.</li>
<li>When transforming data files, we stick to using Picard tools over other tools to avoid subtle incompatibilities.</li>
<li>Finally, when I call default options within a command, follow suit to ensure the same results.    </li>
</ul>
<hr />
<h4>The steps of the workflow are as follows.</h4>
<ol>
<li><a href="#step1">Generate an unmapped BAM (uBAM)</a>
(A) Convert the FASTQ to uBAM and add read group information using FastqToSam
(B1) [Optional] Extract reads in a genomic interval from aligned BAM
(B2) Convert aligned BAM to uBAM and discard problematic records using RevertSam   </li>
<li><a href="#step2">Mark adapter sequences using MarkIlluminaAdapters</a></li>
<li><a href="#step3">Convert uBAM to FASTQ and assign adapter bases low qualities using SamToFastq</a></li>
<li><a href="#step4">Align reads and flag secondary hits using BWA MEM</a></li>
<li>[Optional] <a href="#step5">Pipe steps 3 &amp; 4 and collect alignment metrics</a></li>
<li>[Optional] <a href="#step6">Sort, index and convert alignment to a BAM using SortSam and visualize on IGV</a></li>
<li><a href="#step7">Restore altered data and apply &amp; adjust meta information using MergeBamAlignment</a></li>
</ol>
<hr />
<p><a name="step1"></a></p>
<h3>1. Generate an unmapped BAM (uBAM)</h3>
<p>The goal is to produce an unmapped BAM file with <em>appropriate</em> read group (@RG) information that differentiates not only samples, but also factors that contribute to technical artifacts. To see the read group information for a BAM file, use the following command. </p>
<pre><code class="pre_md">samtools view -H Solexa-272222.bam | grep '@RG'</code class="pre_md"></pre>
<p>This prints the lines starting with @RG within the header. Our tutorial file's single @RG line is shown below. The file has the read group fields required by this workflow as well as extra fields for record keeping. Two read group fields, <code>ID</code> and <code>PU</code>, appropriately differentiate flow cell lane, marked by <code>.2</code>, a factor that contributes to batch effects.  </p>
<pre><code class="pre_md">@RG ID:H0164.2  PL:illumina PU:H0164ALXX140820.2    LB:Solexa-272222    PI:0    DT:2014-08-20T00:00:00-0400 SM:NA12878  CN:BI</code class="pre_md"></pre>
<ul>
<li>GATK's variant discovery workflow requires <code>ID</code>, <code>SM</code> and <code>LB</code> fields and recommends the <code>PL</code> field. </li>
<li>Each <code>@RG</code> line has a unique <code>ID</code> that differentiates read groups. It is the lowest denominator that differentiates factors contributing to technical batch effects and is repeatedly indicated by the <code>RG</code> tag for each read record. Thus, the length of this field contributes to file size. </li>
<li><code>SM</code> indicates sample name and, within a collection of samples, <code>LB</code> indicates if the same sample was sequenced in multiple lanes. See item 8 of <a href="http://gatkforums.broadinstitute.org/discussion/1317/collected-faqs-about-bam-files">Collected FAQs about BAM files</a> for more detail. </li>
<li><code>PU</code> is not required by any GATK tool. If present it is used by BQSR instead of <code>ID</code>. It is required by Picard's AddOrReplaceReadGroups but not FastqToSam. </li>
</ul>
<p>If your sample collection's BAM files lack required fields or do not differentiate pertinent factors within the fields, use Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a> to add or appropriately rename the read group fields.</p>
<p>Here we illustrate how to derive both <code>ID</code> and <code>PU</code> fields from query names. We break down the common portion of two different read query names from the tutorial file. </p>
<pre><code class="pre_md">H0164ALXX140820:2:1101:10003:23460
H0164ALXX140820:2:1101:15118:25288

#Breaking down the common portion of the query names:
H0164____________ # portion of @RG ID and PU fields indicating Illumina flow cell
_____ALXX140820__ # portion of @RG PU field indicating barcode or index in a multiplexed run
_______________:2 # portion of @RG ID and PU fields indicating flow cell lane</code class="pre_md"></pre>
<hr />
<h4>(A)  Convert the FASTQ to uBAM and add read group information using FastqToSam</h4>
<p>Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam">FastqToSam</a> transforms a FASTQ file to unmapped BAM, requires two read group fields and makes optional specification of other read group fields. In the command below we note which fields are required for our workflow. All other read group fields are optional. </p>
<pre><code class="pre_md">java -Xmx8G -jar /seq/software/picard/current/bin/picard.jar FastqToSam \
    FASTQ=snippet_XT_interleaved.fq \ #our single tutorial file contains both reads in a pair 
    OUTPUT=snippet_FastqToSam_PU.bam \
    READ_GROUP_NAME=H0164.2 \ # required; changed from default of A
    SAMPLE_NAME=NA12878 \ # required
    LIBRARY_NAME=Solexa-272222 \ # required 
    PLATFORM_UNIT=H0164ALXX140820.2 \ 
    PLATFORM=illumina \ # recommended
    SEQUENCING_CENTER=BI \ 
    RUN_DATE=2014-08-20T00:00:00-0400</code class="pre_md"></pre>
<p>Some details on select parameters:    </p>
<ul>
<li><code>QUALITY_FORMAT</code> is detected automatically if unspecified.</li>
<li><code>SORT_ORDER</code> by default is queryname.</li>
<li>Specify both <code>FASTQ</code> and <code>FASTQ2</code> for paired reads in separate files. </li>
<li><code>PLATFORM_UNIT</code> is often in run_barcode.lane format. Include if sample is multiplexed.</li>
<li><code>RUN_DATE</code> is in <a href="https://en.wikipedia.org/wiki/ISO_8601">Iso8601 date format</a>.</li>
</ul>
<hr />
<h4>(B1) [Optional] Extract reads in a genomic interval from aligned BAM</h4>
<p>We want to test our reversion process on a subset of the tutorial file before committing to reverting the entire BAM. This process requires the reads in the BAM to be aligned to a reference genome and produces a BAM containing reads from a genomic interval.</p>
<pre><code class="pre_md">java -Xmx8G -jar /path/GenomeAnalysisTK.jar \
    -T PrintReads \ 
    -R /path/human_g1k_v37_decoy.fasta \
    -L 10:90000000-100000000 \ # this is the retained interval
    -I Solexa-272222.bam -o snippet.bam # snippet.bam is newly created</code class="pre_md"></pre>
<ul>
<li>This seems a good time to bring this up. In the command, the <code>-Xmx8G</code> Java option sets the maximum heap size, or memory usage to eight gigabytes. We want to both cap Java's use of memory so the system doesn't slow down as well as allow enough memory for the tool to run without causing an out of memory error. The <code>-Xmx</code> settings we provide here is more than sufficient for most cases. For GATK, 4G is standard, while for Picard less is needed. Some tools, e.g. MarkDuplicates, may require more. I have heard up to16G specified and have also omitted this option for small files. To find a system's default maximum heap size, type <code>java -XX:+PrintFlagsFinal -version</code>, and look for <code>MaxHeapSize</code>. Note that any setting beyond available memory spills to storage and slows a system down. If <a href="https://www.broadinstitute.org/gatk/guide/article?id=1975">multithreading</a>, increase memory proportionately to the number of threads. e.g. if 1G is the minimum required for one thread, then use 2G for two threads.</li>
<li>This step is for our tutorial only. For applying interval lists, e.g. to whole exome data, see <a href="http://gatkforums.broadinstitute.org/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals">When should I use L to pass in a list of intervals</a>.</li>
</ul>
<hr />
<h4>(B2) Convert aligned BAM to uBAM and discard problematic records using RevertSam</h4>
<p>We use Picard's RevertSam to remove alignment information. The resulting unmapped BAM (uBAM) has two uses in this workflow: (1) for processing through the MarkIlluminaAdapters branch of the workflow, and (2) for application of read group, read sequence and other read meta information to the aligned read file in the MergeBamAlignment branch of the workflow. The RevertSam parameters we specify remove information pertaining to previous alignments including program group records and standard alignment flags and tags that would otherwise transfer over in the MergeBamAlignment step. We remove nonstandard alignment tags with the <code>ATTRIBUTE_TO_CLEAR</code> option. For example, we clear the <code>XT</code> tag using this option so that it is free for use by MarkIlluminaAdapters. Our settings also reset <a href="https://broadinstitute.github.io/picard/explain-flags.html">flags</a> to unmapped values, e.g. 77 and 141 for paired reads.  Additionally, we invoke the <code>SANITIZE</code> option to remove reads that cause problems for MarkIlluminaAdapters. Our tutorial's <code>snippet</code> requires such filtering while <code>Solexa-272222</code> does not. </p>
<p>For our particular file, we use the following parameters.</p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar RevertSam \
    I=snippet.bam \
    O=snippet_revert.bam \
    SANITIZE=true \ 
    MAX_DISCARD_FRACTION=0.005 \ # informational; does not affect processing
    ATTRIBUTE_TO_CLEAR=XT \
    ATTRIBUTE_TO_CLEAR=XN \
    ATTRIBUTE_TO_CLEAR=AS \ #Picard release of 9/2015 clears AS by default
    ATTRIBUTE_TO_CLEAR=OC \
    ATTRIBUTE_TO_CLEAR=OP \
    SORT_ORDER=queryname \ #default
    RESTORE_ORIGINAL_QUALITIES=true \ #default
    REMOVE_DUPLICATE_INFORMATION=true \ #default
    REMOVE_ALIGNMENT_INFORMATION=true #default</code class="pre_md"></pre>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee # sets environmental variable for temporary directory</code class="pre_md"></pre>
<p>We change these settings for RevertSam:</p>
<ul>
<li><code>SANITIZE</code> If the BAM file contains problematic reads, such as that might arise from taking a genomic interval of reads (Step 1), then RevertSam's <code>SANTITIZE</code> option removes them. Our workflow's downstream tools will have problems with paired reads with missing mates, duplicated records, and records with mismatches in length of bases and qualities. </li>
<li>
<p><code>MAX_DISCARD_FRACTION</code> is set to a more strict threshold of 0.005 instead of the default 0.01. Whether or not this fraction is reached, the tool informs you of the number and fraction of reads it discards. This parameter asks the tool to additionally inform you of the discarded fraction via an exception as it finishes processing. </p>
<pre><code class="pre_md">Exception in thread "main" picard.PicardException: Discarded 0.947% which is above MAX_DISCARD_FRACTION of 0.500%  </code class="pre_md"></pre>
</li>
<li>
<p><code>ATTRIBUTE_TO_CLEAR</code> is set to clear more than the default standard tags, which are NM, UQ, PG, MD, MQ, SA, MC, and AS tags. The AS tag is removed by default for Picard releases starting 9/2015. Remove all other tags, such as the XT tag needed by MarkIlluminaAdapters, by specifying each with the <code>ATTRIBUTE_TO_CLEAR</code> option. To list all tags within my BAM, I used the command below to get RG, OC, XN, OP, <em>SA</em>, <em>MD</em>, <em>NM</em>, <em>PG</em>, <em>UQ</em>, <em>MC</em>, <em>MQ</em>, <em>AS</em>, XT, and <em>OQ</em> tags. Those removed by default and by <code>RESTORE_ORIGINAL_QUALITIES</code> are italicized. See your aligner's documentation and the <a href="http://samtools.sourceforge.net/SAM1.pdf">Sequence Alignment/Map Format Specification</a> for descriptions of tags.   </p>
<pre><code class="pre_md">samtools view input.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}' </code class="pre_md"></pre>
</li>
</ul>
<p>Some comments on options kept at default:</p>
<ul>
<li><code>SORT_ORDER</code>=queryname
For paired read files, because each read in a pair has the same query name, sorting results in interleaved reads. This means that reads in a pair are listed consecutively within the same file. We make sure to alter the previous sort order. Coordinate sorted reads result in the aligner incorrectly estimating insert size from blocks of paired reads as they are not randomly distributed. </li>
<li><code>RESTORE_ORIGINAL_QUALITIES</code>=true
Restoring original base qualities to the QUAL field requires OQ tags listing original qualities. The OQ tag uses the same encoding as the QUAL field, e.g. ASCII Phred-scaled base quality+33 for tutorial data. After restoring the QUAL field, RevertSam removes the tag.</li>
<li><code>REMOVE_ALIGNMENT_INFORMATION</code>=true will remove program group records and alignment information. It also invokes the default <code>ATTRIBUTE_TO_CLEAR</code> parameter which removes standard alignment tags.</li>
</ul>
<p>For snippet.bam, <code>SANITIZE</code> removes 25,909 out of 2,735,539 (0.947%) reads, leaving us with 2,709,630 reads. The intact BAM retains all reads. The example shows a read pair before and after RevertSam. </p>
<pre><code class="pre_md">#original BAM
H0164ALXX140820:2:1101:10003:23460  83  10  91515318    60  151M    =   91515130    -339    CCCATCCCCTTCCCCTTCCCTTTCCCTTTCCCTTTTCTTTCCTCTTTTAAAGAGACAAGGTCTTGTTCTGTCACCCAGGCTGGAATGCAGTGGTGCAGTCATGGCTCACTGCCGCCTCAGACTTCAGGGCAAAAGCAATCTTTCCAGCTCA :&lt;&lt;=&gt;@AAB@AA@AA&gt;6@@A:&gt;,*@A@&lt;@??@8?9&gt;@==8?:?@?;?:&gt;&lt;??@&gt;==9?&gt;8&gt;@:?&gt;&gt;=&gt;;&lt;==&gt;&gt;;&gt;?=?&gt;&gt;=&lt;==&gt;&gt;=&gt;9&lt;=&gt;??&gt;?&gt;;8&gt;?&gt;&lt;?&lt;=:&gt;&gt;&gt;;4&gt;=&gt;7=6&gt;=&gt;&gt;=&gt;&lt;;=;&gt;===?=&gt;=&gt;&gt;?9&gt;&gt;&gt;&gt;??==== MC:Z:60M91S MD:Z:151    PG:Z:MarkDuplicates RG:Z:H0164.2    NM:i:0  MQ:i:0  OQ:Z:&lt;FJFFJJJJFJJJJJF7JJJ&lt;F--JJJFJJJJ&lt;J&lt;FJFF&lt;JAJJJAJAJFFJJJFJAFJAJJAJJJJJFJJJJJFJJFJJJJFJFJJJJFFJJJJJJJFAJJJFJFJFJJJFFJJJ&lt;J7JJJJFJ&lt;AFAJJJJJFJJJJJAJFJJAFFFFA    UQ:i:0  AS:i:151
H0164ALXX140820:2:1101:10003:23460  163 10  91515130    0   60M91S  =   91515318    339 TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCCCCTCTCCCACCCCTCTCTCCCCCCCTCCCACCC :0;.=;8?7==?794&lt;&lt;;:&gt;769=,&lt;;0:=&lt;0=:9===/,:-==29&gt;;,5,98=599;&lt;=########################################################################################### SA:Z:2,33141573,-,37S69M45S,0,1;    MC:Z:151M   MD:Z:48T4T6 PG:Z:MarkDuplicates RG:Z:H0164.2    NM:i:2  MQ:i:60 OQ:Z:&lt;-&lt;-FA&lt;F&lt;FJF&lt;A7AFAAJ&lt;&lt;AA-FF-AJF-FA&lt;AFF--A-FA7AJA-7-A&lt;F7&lt;&lt;AFF###########################################################################################    UQ:i:49 AS:i:50

#after RevertSam (step 1.B2)
H0164ALXX140820:2:1101:10003:23460  77  *   0   0   *   *   0   0   TGAGCTGGAAAGATTGCTTTTGCCCTGAAGTCTGAGGCGGCAGTGAGCCATGACTGCACCACTGCATTCCAGCCTGGGTGACAGAACAAGACCTTGTCTCTTTAAAAGAGGAAAGAAAAGGGAAAGGGAAAGGGAAGGGGAAGGGGATGGG AFFFFAJJFJAJJJJJFJJJJJAFA&lt;JFJJJJ7J&lt;JJJFFJJJFJFJFJJJAFJJJJJJJFFJJJJFJFJJJJFJJFJJJJJFJJJJJAJJAJFAJFJJJFFJAJAJJJAJ&lt;FFJF&lt;J&lt;JJJJFJJJ--F&lt;JJJ7FJJJJJFJJJJFFJF&lt; RG:Z:H0164.2
H0164ALXX140820:2:1101:10003:23460  141 *   0   0   *   *   0   0   TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCCCCTCTCCCACCCCTCTCTCCCCCCCTCCCACCC &lt;-&lt;-FA&lt;F&lt;FJF&lt;A7AFAAJ&lt;&lt;AA-FF-AJF-FA&lt;AFF--A-FA7AJA-7-A&lt;F7&lt;&lt;AFF########################################################################################### RG:Z:H0164.2</code class="pre_md"></pre>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step2"></a></p>
<h3>2. Mark adapter sequences using MarkIlluminaAdapters</h3>
<p>Previously we cleared the XT tag from our BAM so Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters">MarkIlluminaAdapters</a> can use it to mark adapter sequences. SamToFastq (step 4) will use these in turn to assign low base quality scores to the adapter bases, effectively removing their contribution to read alignment and alignment scoring metrics. For the tutorial data, adapter sequences have already been removed from the <em>beginning</em> of reads. We want to additionally effectively remove any adapter sequences at the <em>ends</em> of reads arising from read-through to adapters in read pairs with shorter inserts. </p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar MarkIlluminaAdapters \
    I=snippet_revert.bam \
    O=snippet_revertmark.bam \
    M=snippet_revertmark.metrics.txt \ #naming required
    TMP_DIR=/path/shlee # optional to process large files</code class="pre_md"></pre>
<ul>
<li>By default, the tool uses Illumina adapter sequences. This is sufficient for our tutorial data. Specify other adapter sequences as outlined in the <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters">tool documentation</a>.</li>
<li>Only reads with adapter sequence are marked with the tag in XT:i:[#] format, where # denotes the starting position of the adapter sequence.  </li>
</ul>
<p>The example shows a read pair marked with the XT tag by MarkIlluminaAdapters. This is a different pair than shown previously as <code>H0164ALXX140820:2:1101:10003:23460</code> reads do not contain adapter sequence. The insert region sequences for the reads overlap by a length corresponding approximately to the XT tag value. The same read pair is shown after SamToFastq transformation, where adapter sequence base quality scores have been set to 2 (# symbol), and after MergeBamAlignment, which restores original base quality scores. </p>
<pre><code class="pre_md">#after MarkIlluminaAdapters (step 2)
H0164ALXX140820:2:1101:15118:25288  77  *   0   0   *   *   0   0   
ACCTGCCTCAGCCTCCCAAAGTGCTGGGATTATAGGTATGTGTCACCACACCCAGCCAAGTATACTCACATTGTCGTGCAACCAAACTCCAGAACTTTTTCATCTTAAAGAATCAAGGTTTTTTATTGTTTACTTTATTACTTATTTATTT 
AFFFFFJJFJFAAJJFFJJFJFJ&lt;FJJJJJJF&lt;JJJFFJJAF7JJJAAF7AJJFJFJFFJ--A-FAJA-F&lt;J7A--AFJ7AJ7AJ-FJ7-JJJ-F-J---7J---7FF-JAJJ&lt;A7JFAFAA7--FF----AF-7&lt;JF&lt;JFA-7&lt;F-FF-J RG:Z:H0164.2    XT:i:63
H0164ALXX140820:2:1101:15118:25288  141 *   0   0   *   *   0   0   
GTCATGGCTGGACGCAGTGGCTCATACCTGTAATCCCAGCACTTTTGGAGGCTGAGGCAGGTAGATCGGAAGCGCCTCGTGTAGGGAGAGAGGGTTAACAAAAATGTAGATACCGGAGGTCGCCGTAAAATAAAAAAGTAGCAAGGAGTAG 
AAFFFJJJJJAJJJJJFJJJJ&lt;JFJJJJJJJJFJJJJFJ&lt;FJJJJAJJJJJJJJFJJJ7JJ--JJJ&lt;J&lt;-FJ7F--&lt;-J7--7AJJA-J------J7F&lt;-77--F--FFJ---J-J-J--A-7&lt;&lt;----J-7-J-FJ--J--FA####### RG:Z:H0164.2    XT:i:63

#after SamToFastq (step 3)
@H0164ALXX140820:2:1101:15118:25288/1
ACCTGCCTCAGCCTCCCAAAGTGCTGGGATTATAGGTATGTGTCACCACACCCAGCCAAGTATACTCACATTGTCGTGCAACCAAACTCCAGAACTTTTTCATCTTAAAGAATCAAGGTTTTTTATTGTTTACTTTATTACTTATTTATTT
+
AFFFFFJJFJFAAJJFFJJFJFJ&lt;FJJJJJJF&lt;JJJFFJJAF7JJJAAF7AJJFJFJFFJ--#########################################################################################
@H0164ALXX140820:2:1101:15118:25288/2
GTCATGGCTGGACGCAGTGGCTCATACCTGTAATCCCAGCACTTTTGGAGGCTGAGGCAGGTAGATCGGAAGCGCCTCGTGTAGGGAGAGAGGGTTAACAAAAATGTAGATACCGGAGGTCGCCGTAAAATAAAAAAGTAGCAAGGAGTAG
+
AAFFFJJJJJAJJJJJFJJJJ&lt;JFJJJJJJJJFJJJJFJ&lt;FJJJJAJJJJJJJJFJJJ7JJ-#########################################################################################

#after MergeBamAlignment (step 7)
H0164ALXX140820:2:1101:15118:25288  99  10  99151971    60  151M    =   99152350    440 
ACCTGCCTCAGCCTCCCAAAGTGCTGGGATTATAGGTATGTGTCACCACACCCAGCCAAGTATACTCACATTGTCGTGCAACCAAACTCCAGAACTTTTTCATCTTAAAGAATCAAGGTTTTTTATTGTTTACTTTATTACTTATTTATTT
AFFFFFJJFJFAAJJFFJJFJFJ&lt;FJJJJJJF&lt;JJJFFJJAF7JJJAAF7AJJFJFJFFJ--A-FAJA-F&lt;J7A--AFJ7AJ7AJ-FJ7-JJJ-F-J---7J---7FF-JAJJ&lt;A7JFAFAA7--FF----AF-7&lt;JF&lt;JFA-7&lt;F-FF-J MC:Z:90S61M MD:Z:74T10T3A37T23  PG:Z:bwamem RG:Z:H0164.2    NM:i:4  MQ:i:60 UQ:i:48 AS:i:131    XS:i:40
H0164ALXX140820:2:1101:15118:25288  147 10  99152350    60  90S61M  =   99151971    -440
CTACTCCTTGCTACTTTTTTATTTTACGGCGACCTCCGGTATCTACATTTTTGTTAACCCTCTCTCCCTACACGAGGCGCTTCCGATCTACCTGCCTCAGCCTCCAAAAGTGCTGGGATTACAGGTATGAGCCACTGCGTCCAGCCATGAC 
#######AF--J--JF-J-7-J----&lt;&lt;7-A--J-J-J---JFF--F--77-&lt;F7J------J-AJJA7--7J-&lt;--F7JF-&lt;J&lt;JJJ--JJ7JJJFJJJJJJJJAJJJJF&lt;JFJJJJFJJJJJJJJFJ&lt;JJJJFJJJJJAJJJJJFFFAA MC:Z:151M   MD:Z:61 PG:Z:bwamem RG:Z:H0164.2    NM:i:0  MQ:i:60 UQ:i:0  AS:i:61 XS:i:50</code class="pre_md"></pre>
<p>Snippet_revertmark.bam marks 5,810 reads (0.21%) with XT, while Solexa-272222_revertmark.bam marks 3,236,552 reads (0.39%). We plot the metrics data using RStudio.
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/6e/2da8652875645713c23e45fddef790.png" height="230" border="9"/> <img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a1/59e1837963fe37d0577d466f3c56b2.png"height="230" border="9" /></p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3"></a></p>
<h3>3. Convert BAM to FASTQ using SamToFastq</h3>
<p>Picard's SamToFastq takes read identifiers, read sequences, and base quality scores to write a Sanger FASTQ format file. We use additional options to effectively remove adapter sequences previously marked with the XT tag. All extant meta data, i.e. alignment information, flags and tags, are purged in this transformation. </p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar SamToFastq \
    I=snippet_revertmark.bam \
    FASTQ=snippet_XT_interleaved.fq \
    CLIPPING_ATTRIBUTE=XT \
    CLIPPING_ACTION=2 \
    INTERLEAVE=true \ 
    NON_PF=true \
    TMP_DIR=/path/shlee # optional to process large files         </code class="pre_md"></pre>
<ul>
<li>
<p>By specifying <code>CLIPPING_ATTRIBUTE</code>=XT and <code>CLIPPING_ACTION</code>=2, SamToFastq changes the quality scores of bases marked by XT to two--a rather low score in the Phred scale. This effectively removes the adapter portion of sequences from contributing to read alignment and alignment scoring metrics. This reassignment is temporary as we will restore the original base quality scores after alignment in step 7.</p>
</li>
<li>
<p>For our paired reads sample we set SamToFastq's <code>INTERLEAVE</code> to true. During the conversion to FASTQ format, the query name of the reads in a pair are marked with /1 or /2 and paired reads are retained in the same FASTQ file.</p>
<p><a href="http://bio-bwa.sourceforge.net/bwa.shtml">BWA aligner</a> accepts interleaved FASTQ files given the <code>-p</code> option. This command indicates that the i-th and the (i+1)-th reads constitute a read pair. </p>
</li>
<li>We change the <code>NON_PF</code>, aka <code>INCLUDE_NON_PF_READS</code>, option from default to true. SamToFastq will then retain reads marked by what <a href="https://github.com/samtools/hts-specs/issues/85">some consider archaic 0x200 flag bit</a> that denotes reads that do not pass quality controls. These reads are also known as failing platform or vendor quality checks. Our tutorial data do not contain such reads and we call out this option for illustration only.</li>
</ul>
<h4>[Optional] Compress the FASTQ using gzip</h4>
<p>This step is optional. The step is irrelevant if you pipe steps 3 and 4, as we outline in step 5.  </p>
<p>BWA handles both FASTQ and gzipped FASTQ files natively--that is, BWA works on both file types directly. Thus, this step is optional. Compress the FASTQ file using the UNIX gzip utility. </p>
<pre><code class="pre_md">gzip snippet_XT_interleaved.fq #replaces the file with snippet_XT_interleaved.fq.gz</code class="pre_md"></pre>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step4"></a></p>
<h3>4. Align reads and flag secondary hits using BWA MEM</h3>
<p>GATK's variant discovery workflow recommends Burrows-Wheeler Aligner's maximal exact matches (BWA MEM) algorithm (<a href="http://arxiv.org/abs/1303.3997">Li 2013 reference</a>; <a href="http://bioinformatics.oxfordjournals.org/content/30/20/2843.long">Li 2014 benchmarks</a>; <a href="http://bio-bwa.sourceforge.net/">homepage</a>; <a href="http://bio-bwa.sourceforge.net/bwa.shtml">manual</a>). BWA MEM is suitable for aligning high-quality long reads ranging from 70 bp to 1 Mbp against a large reference genome such as the human genome.  </p>
<ul>
<li>We use BWA v 0.7.7.r441, the same aligner used by the Broad's Genomics Platform as of this writing (9/2015).</li>
<li>Alignment is a compute intensive process. For faster processing, use a reference genome with decoy sequences, also called a <a href="http://www.cureffi.org/2013/02/01/the-decoy-genome/">decoy genome</a>. For example, the Broad's Genomics Platform uses an Hg19/GRCh37 reference sequence that includes Ebstein-Barr virus (EBV) sequence to soak up reads that fail to align to the human reference that the aligner would otherwise spend an inordinate amount of time trying to align as split reads. <a href="https://www.broadinstitute.org/gatk/guide/article.php?id=1213">GATK's resource bundle</a> provides a standard decoy genome from the <a href="http://www.1000genomes.org/">1000 Genomes Project</a>.</li>
<li>Aligning our <code>snippet</code> reads from a genomic interval against either a portion or the whole genome is not equivalent to aligning our entire file and taking a new <code>slice</code> from the same genomic interval. </li>
</ul>
<p><strong>Index the reference genome file for BWA.</strong> Indexing is specific to algorithms. To index the human genome for BWA, we apply BWA's <code>index</code> function on the reference genome file, e.g. <code>human_g1k_v37_decoy.fasta</code>. This produces five index files with the extensions <code>amb</code>, <code>ann</code>, <code>bwt</code>, <code>pac</code> and <code>sa</code>. </p>
<pre><code class="pre_md">bwa index -a bwtsw human_g1k_v37_decoy.fasta</code class="pre_md"></pre>
<p><strong>Align using BWA MEM.</strong> The tool automatically locates the index files within the same folder as the reference FASTA file. In the alignment command, <code>&gt;</code> denotes the aligned file. </p>
<ul>
<li>The aligned file is in SAM format even if given a BAM extension and retains the sort order of the FASTQ file. Thus, our aligned tutorial file remains sorted by query name. </li>
<li>
<p>BWA automatically creates a program group record (@PG) in the header that gives the ID, group name, group version, and command line information. </p>
<p>/path/bwa mem -M -t 7 -p \
/path/Homo_sapiens_assembly19.fasta \ #reference genome
Solexa-272222_interleavedXT.fq &gt; Solexa-272222_markXT_aln.sam </p>
</li>
</ul>
<p>We invoke three options in the command. </p>
<ul>
<li><code>-M</code> to flag shorter split hits as secondary.
This is optional for Picard compatibility. However, if we want MergeBamAlignment to reassign proper pair alignments, we need to mark secondary alignments.  </li>
<li><code>-p</code> to indicate the given file contains interleaved paired reads.</li>
<li>
<p><code>-t</code> followed by a number for the <em>additional</em> number of processor threads to use concurrently. Check your server or system's total number of threads with the following command.</p>
<pre><code class="pre_md">getconf _NPROCESSORS_ONLN #thanks Kate</code class="pre_md"></pre>
</li>
</ul>
<p>MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, the point of this workflow is to take advantage of the features offered by MergeBamAlignment that allow for the scalable, <em>lossless</em> operating procedure practiced by Broad's Genomics Platform and to produce comparable metrics.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step5"></a></p>
<h3>5. [Optional] Pipe steps 3 &amp; 4 and collect alignment metrics</h3>
<p><strong>Piping processes saves time and space.</strong> Our tutorial's resulting SAM file is small enough to easily view, manipulate and store. For larger data, however, consider using <a href="https://en.wikipedia.org/wiki/Pipeline_(Unix)">Unix pipelines</a>. Piping allows streaming data in the processor's input-output (I/O) device directly to the next process for efficient processing and storage. We recommend piping steps 3 and 4 so as to avoid rereading and storing the large intermediate FASTQ file. </p>
<p>You may additionally extend piping to include step 6's SortSam. Steps 3-4-6 are piped in the example command below to generate an aligned BAM file and index. [For the larger file, I couldn't pipe Step 7's MergeBamAlignment.]</p>
<pre><code class="pre_md">#overview of command structure
[step 3's SamToFastq] | [step 4's bwa mem] | [step 6's SortSam]

#for our file  
java -Xmx8G -jar /path/picard.jar SamToFastq I=snippet_revertmark.bam \
    FASTQ=/dev/stdout \
    CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
    TMP_DIR=/path/shlee | \ 
    /path/bwa mem -M -t 7 -p /path/Homo_sapiens_assembly19.fasta \
    /dev/stdin | \  #to stop piping here, add '&gt; snippet_piped.sam'
    java -Xmx8G -jar /path/picard.jar SortSam \
    INPUT=/dev/stdin \
    OUTPUT=snippet_piped.bam \
    SORT_ORDER=coordinate CREATE_INDEX=true \
    TMP_DIR=/path/shlee</code class="pre_md"></pre>
<p><strong>Calculate alignment metrics using Picard tools.</strong> Picard offers a variety of metrics collecting tools, e.g. <a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics">CollectAlignmentSummaryMetrics</a>, <a href="http://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics">CollectWgsMetrics</a> and <a href="http://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics">CollectInsertSizeMetrics</a>. Some tools give more detailed metrics if given the reference sequence. See <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html">Picard for metrics definitions</a>. Metrics calculations will differ if run on the BAM directly from alignment (BWA) versus on the merged BAM (MergeBamAlignment). See [link--get from G] for guidelines on when to run tools.  </p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar CollectAlignmentSummaryMetrics \
    R=/path/Homo_sapiens_assembly19.fasta \
    INPUT=slice.bam \
    OUTPUT=slice_bam_metrics.txt \
    TMP_DIR=/path/shlee # optional to process large files</code class="pre_md"></pre>
<p>For example, percent chimeras is a calculated metric. Our tutorial alignment of the whole data set gives 0.019% (BWA) or 0.0034% (MergeBamAlignment) chimeric paired reads. The genomic interval defined in step 1 reports 0.0032% chimeric paired reads. In contrast, the aligned <code>snippet</code> gives 0.0012% (BWA) or 0.00002% (MergeBamAlignment) chimeric paired reads. This illustrates in part the differences I alluded to at the beginning of step 4.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step6"></a></p>
<h3>6. [Optional] Sort, index and convert alignment to a BAM using SortSam and visualize on IGV</h3>
<p><strong>Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#SortSam">SortSam</a> sorts, indexes and converts between SAM and BAM formats.</strong> For file manipulations and to view aligned reads using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>, the SAM or BAM file must be coordinate-sorted and indexed. Some Picard tools, such as MergeBamAlignment in step 7, by default coordinate sort and can use the standard <code>CREATE_INDEX</code> option. If you didn't create an index in step 7, or want to convert to BAM and index the alignment file from step 4, then use Picard's SortSam. The index file will have an <code>sai</code> or <code>bai</code> extension depending on the specified format.</p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar SortSam \
    INPUT=Solexa-272222_markXT_aln.sam \ 
    OUTPUT=Solexa-272222_markXT_aln.bam \ #extension here specifies format conversion
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \ # a standard option for Picard commands
    TMP_DIR=/path/shlee # optional to process large files</code class="pre_md"></pre>
<p><strong>View aligned reads using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>.</strong> Of the multiple IGV versions, the Java Web Start <code>jnlp</code> version allows the highest memory, as of this writing 10 GB for machines with 64-bit Java. </p>
<ul>
<li>To run the <code>jnlp</code> version of IGV, you may need to adjust your system's <em>Java Control Panel</em> settings, e.g. enable Java content in the browser. Also, when first opening the <code>jnlp</code>, overcome Mac OS X's gatekeeper function by right-clicking the saved <code>jnlp</code> and selecting <em>Open</em> <em>with Java Web Start</em>. </li>
<li>Load the appropriate reference genome. For our tutorial this is <em>Human (b37)</em>. </li>
<li>Go to <em>View</em>&gt;<em>Preferences</em> and make sure the settings under the <em>Alignments</em> tab allows you to view reads of interest, e.g. duplicate reads. Default settings are tuned to genomic sequence libraries. Right-click on a track to access a menu of additional viewing modes. See <a href="http://www.broadinstitute.org/igv/AlignmentData">Viewing Alignments</a> in IGV  documentation for details.</li>
<li>Go to <em>File</em>&gt;<em>Load from</em> and either load alignments from <em>File</em> or <em>URL</em>. </li>
</ul>
<p>Here, IGV displays our example chimeric pair, <code>H0164ALXX140820:2:1101:10003:23460</code> at its alignment loci. BWA's secondary alignment designation causes the mates on chromosome 10 to display as unpaired in IGV's paired view. MergeBamAlignment corrects for this when it switches the secondary alignment designation. Mates display as paired on chromosome 10. </p>
<p>Visualizing alignments in such a manner makes apparent certain convergent information. For example, we see that the chimeric region on chromosome 2 is a low complexity GC-rich region, apparent by the predominantly yellow coloring (representing guanine) of the reference region. We know there are many multimapping reads because reads with MAPQ score of zero are filled in white versus gray, and the region is down-sampled, as indicated by the underscoring in the log-scaled coverage chart. We can infer reads in this chromosome 2 region are poorly mapped based on the region's low complexity, depth of reads and prevalence of low MAPQ reads. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a4/d4b30ab2ad0b6595539ee25e115ffd.png"  border="9"/>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step7"></a></p>
<h3>7. Restore altered data and apply &amp; adjust meta information using MergeBamAlignment</h3>
<p>Our alignment file lacks read group information and certain tags, such as the mate CIGAR (MC) tag. It has hard-clipped sequences and altered base qualities. The alignment also has some mapping artifacts we would like to correct for accounting congruency. Finally, the alignment records require coordinate sorting and indexing. </p>
<p>We use Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a> to address all of these needs to produce a <em>raw</em> BAM file that is ready for GATK's variant discovery workflow. MergeBamAlignment takes metadata from a SAM or BAM file of unmapped reads (uBAM) and merges it with a SAM or BAM file containing alignment records for a <em>subset</em> of those reads. Metadata include read group information, read sequences, base quality scores and tags. The tool applies read group information from the uBAM and retains the program group information from the aligned file. In restoring original sequences, MergeBamAlignment adjusts CIGAR strings from hard-clipped to soft-clipped. The tool adjusts <a href="https://broadinstitute.github.io/picard/explain-flags.html">flag</a> values, e.g. changes primary alignment designations according to a user-specified strategy, for desired congruency. Optional parameters allow introduction of additional metadata, e.g. user-specified program group information or nonstandard aligner-generated tags. If the alignment file is missing reads present in the unaligned file, these are retained as unaligned records. Finally, alignment records are coordinate sorted, meaning they are ordered by chromosomal mapping position.</p>
<ul>
<li>To simply edit read group information, see Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a>. To simply concatenate read records into one file, use Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#GatherBamFiles">GatherBamFiles</a>. An advantage of using MergeBamAlignment over AddOrReplaceReadGroups is the ability to transfer mixed read groups to reads in a single file.</li>
<li>Consider what <code>PRIMARY_ALIGNMENT_STRATEGY</code> option best suits your samples. MergeBamAlignment applies this strategy to a read for which the aligner has provided more than one primary alignment, and for which one is designated primary by virtue of another record being marked secondary. MergeBamAlignment considers and switches only existing primary and secondary designations. </li>
<li>MergeBamAlignment retains secondary alignments with the <code>INCLUDE_SECONDARY_ALIGNMENTS</code> parameter. It may be that alignments marked as secondary are truer to biology or at least reveal useful insight.</li>
</ul>
<p>A read with multiple alignment records may map to multiple loci or may be chimeric--that is, splits the alignment. It is possible for an aligner to produce multiple alignments as well as multiple primary alignments, e.g. in the case of a linear alignment set of split reads. When one alignment, or alignment set in the case of chimeric read records, is designated primary, others are designated either secondary or supplementary. Invoking the <code>-M</code> option, we had BWA mark the record with the longest aligning section of split reads as primary and all other records as secondary. MergeBamAlignment further adjusts this secondary designation and other flags, e.g. read mapped in proper pair and mate unmapped flags, to fix mapping artifacts. We only note some changes made by MergeBamAlignment to our tutorial data and by no means comprehensively list its features.</p>
<pre><code class="pre_md">java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
    R=/path/Homo_sapiens_assembly19.fasta \ 
    UNMAPPED_BAM=Solexa-272222_revertclean.bam \ 
    ALIGNED_BAM=Solexa-272222_markXT_aln.sam \
    O=Solexa-272222_merge_IGV_raw.bam \ #output file name in SAM or BAM format
    CREATE_INDEX=true \ #standard option for any Picard command
    ADD_MATE_CIGAR=true \ #default; adds MC tag
    CLIP_ADAPTERS=false \ #changed from default
    CLIP_OVERLAPPING_READS=true \ #default; soft-clips ends so mates do not overlap
    INCLUDE_SECONDARY_ALIGNMENTS=true \ #default
    MAX_INSERTIONS_OR_DELETIONS=-1 \ #changed to allow any number of insertions or deletions
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \ #changed from default BestMapq
    ATTRIBUTES_TO_RETAIN=XS \ #specify multiple times to retain alignment tags starting with X, Y, or Z 
    TMP_DIR=/path/shlee #optional to process large files</code class="pre_md"></pre>
<p>You need not invoke <code>PROGRAM</code> options as BWA's program group information is sufficient and transfer from the alignment during the merging. If, for whatever reason, you need to apply program group information by a different means, then use MergeBamAlignment to assign each of the following program group options. Example information is given. </p>
<pre><code class="pre_md">    PROGRAM_RECORD_ID=bwa \
    PROGRAM_GROUP_NAME=bwamem \
    PROGRAM_GROUP_VERSION=0.7.7-r441 \
    PROGRAM_GROUP_COMMAND_LINE='/path/bwa mem -M -t 7 -p /path/Homo_sapiens_assembly19.fasta Solexa-272222_interleavedXT.fq &gt; Solexa-272222_markXT_aln.sam' \ </code class="pre_md"></pre>
<p>In the command, we change <code>CLIP_ADAPTERS</code>, <code>MAX_INSERTIONS_OR_DELETIONS</code> and <code>PRIMARY_ALIGNMENT_STRATEGY</code> values from default, and invoke other optional parameters.</p>
<ul>
<li>The path to the reference FASTA given by <code>R</code> should also contain the corresponding sequence dictionary with the same base name and extension <code>.dict</code>. Create a sequence dictionary using Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary">CreateSequenceDictionary</a>.</li>
<li><code>CLIP_ADAPTERS</code>=false leaves reads unclipped.</li>
<li><code>MAX_INSERTIONS_OR_DELETIONS</code>=-1 retains reads irregardless of the number of insertions and deletions. Default is 1.</li>
<li><code>PRIMARY_ALIGNMENT_STRATEGY</code>=MostDistant marks primary alignments based on the alignment <em>pair</em> with the largest insert size. This strategy is based on the premise that of chimeric sections of a read aligning to consecutive regions, the alignment giving the largest insert size with the mate gives the most information.</li>
<li><code>ATTRIBUTES_TO_RETAIN</code> is specified to carryover the XS tag from the alignment, which for BWA reports suboptimal alignment scores. The XS tag in not necessary for our workflow. We retain it to illustrate that the tool only carries over select alignment information unless specified otherwise. For our tutorial data, this is the only additional unaccounted tag from the alignment. [IDK if this tag is used downstream; need to confirm I can keep this.]</li>
<li>Because we have left the <code>ALIGNER_PROPER_PAIR_FLAGS</code> parameter at the default false value, MergeBamAlignment may reassign <em>proper pair</em> designations made by the aligner. </li>
<li>By default the merged file is coordinate sorted. We set <code>CREATE_INDEX</code>=true to additionally create the <code>bai</code> index.</li>
</ul>
<p>Original base quality score restoration is illustrated in Step 3. The following example shows a read pair for which MergeBamAlignment adjusts multiple other information fields. The query name is listed thrice because we have paired reads where one of the reads has two alignment loci, on chromosome 2 and on chromosome 10. The mate is mapped with high MAPQ to chromosome 10. The two loci align 69 and 60 nucleotide regions, respectively, and the aligned regions coincide by 15 bases. A good portion of the chromosome 2 aligned region has low base quality scores. The <code>NM</code> tag indicates that the chromosome 2 alignment requires one change to match the reference, while the chromosome 10 read requires two changes and this is also reflected in the <code>MD</code> tags that provide the mismatching positions. When tallying alignment scores, given by the <code>AS</code> tag, aligners penalize mismatching positions, here apparently by five points per mismatch, e.g. 60 matches minus two mismatches multiplied by five gives an alignment score of 50. Both read records have values for the <code>XS</code> (suboptimal alignment score) and <code>SA</code> (chimeric alignment) tags that indicate a split read. Flag values, set by BWA, indicate the chromosome 2 record is primary and the chromosome 10 record is secondary. </p>
<pre><code class="pre_md">#aligned reads from step 4
H0164ALXX140820:2:1101:10003:23460  177 2   33141435    0   37S69M45S   10  91515318    0   
GGGTGGGAGGGGGGGAGAGAGGGGTGGGAGAGGGGAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGAAGGAAAGGAGGGAGGGAGGGAGCAAGGAAGGAAGGAAGGAAAGA ###########################################################################################FFA&lt;&lt;7F&lt;A-7-AJA7AF-A--FFA&lt;AF-FJA-FF-AA&lt;&lt;JAAFA7A&lt;FJF&lt;F&lt;AF-&lt;-&lt; NM:i:1  MD:Z:51G17  AS:i:64 XS:i:64 SA:Z:10,91515130,+,60M91S,0,2;

H0164ALXX140820:2:1101:10003:23460  417 10  91515130    0   60M91H  =   91515318    339 
TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCC    &lt;-&lt;-FA&lt;F&lt;FJF&lt;A7AFAAJ&lt;&lt;AA-FF-AJF-FA&lt;AFF--A-FA7AJA-7-A&lt;F7&lt;&lt;AFF    NM:i:2  MD:Z:48T4T6 AS:i:50 XS:i:36 SA:Z:2,33141435,-,37S69M45S,0,1;

H0164ALXX140820:2:1101:10003:23460  113 10  91515318    60  151M    2   33141435    0
CCCATCCCCTTCCCCTTCCCTTTCCCTTTCCCTTTTCTTTCCTCTTTTAAAGAGACAAGGTCTTGTTCTGTCACCCAGGCTGGAATGCAGTGGTGCAGTCATGGCTCACTGCCGCCTCAGACTTCAGGGCAAAAGCAATCTTTCCAGCTCA &lt;FJFFJJJJFJJJJJF7JJJ&lt;F--JJJFJJJJ&lt;J&lt;FJFF&lt;JAJJJAJAJFFJJJFJAFJAJJAJJJJJFJJJJJFJJFJJJJFJFJJJJFFJJJJJJJFAJJJFJFJFJJJFFJJJ&lt;J7JJJJFJ&lt;AFAJJJJJFJJJJJAJFJJAFFFFA NM:i:0  MD:Z:151    AS:i:151    XS:i:0

#after merging (step 7)
H0164ALXX140820:2:1101:10003:23460  409 2   33141435    0   37S69M45S   =   33141435    0   
GGGTGGGAGGGGGGGAGAGAGGGGTGGGAGAGGGGAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGAAGGAAAGGAGGGAGGGAGGGAGCAAGGAAGGAAGGAAGGAAAGA ###########################################################################################FFA&lt;&lt;7F&lt;A-7-AJA7AF-A--FFA&lt;AF-FJA-FF-AA&lt;&lt;JAAFA7A&lt;FJF&lt;F&lt;AF-&lt;-&lt; SA:Z:10,91515130,+,60M91S,0,2;  MD:Z:51G17  PG:Z:bwamem RG:Z:H0164.2    NM:i:1  UQ:i:2  AS:i:64 XS:i:64

H0164ALXX140820:2:1101:10003:23460  163 10  91515130    0   60M91S  =   91515318    339 
TCTTTCCTTCCTTCCTTCCTTGCTCCCTCCCTCCCTCCTTTCCTTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCCCCTCTCCCACCCCTCTCTCCCCCCCTCCCACCC ###########################################################################################FFA&lt;&lt;7F&lt;A-7-AJA7AF-A--FFA&lt;AF-FJA-FF-AA&lt;&lt;JAAFA7A&lt;FJF&lt;F&lt;AF-&lt;-&lt; SA:Z:2,33141435,-,37S69M45S,0,1;    MC:Z:151M   MD:Z:48T4T6 PG:Z:bwamem RG:Z:H0164.2    NM:i:2  MQ:i:60 UQ:i:4  AS:i:50 XS:i:36

H0164ALXX140820:2:1101:10003:23460  83  10  91515318    60  151M    =   91515130    -339    
CCCATCCCCTTCCCCTTCCCTTTCCCTTTCCCTTTTCTTTCCTCTTTTAAAGAGACAAGGTCTTGTTCTGTCACCCAGGCTGGAATGCAGTGGTGCAGTCATGGCTCACTGCCGCCTCAGACTTCAGGGCAAAAGCAATCTTTCCAGCTCA &lt;FJFFJJJJFJJJJJF7JJJ&lt;F--JJJFJJJJ&lt;J&lt;FJFF&lt;JAJJJAJAJFFJJJFJAFJAJJAJJJJJFJJJJJFJJFJJJJFJFJJJJFFJJJJJJJFAJJJFJFJFJJJFFJJJ&lt;J7JJJJFJ&lt;AFAJJJJJFJJJJJAJFJJAFFFFA MC:Z:60M91S MD:Z:151    PG:Z:bwamem RG:Z:H0164.2    NM:i:0  MQ:i:0  UQ:i:0  AS:i:151    XS:i:0</code class="pre_md"></pre>
<ul>
<li>For the read with two alignments, the aligner hard-clipped the alignment on chromosome 10 giving a CIGAR string of 60M91H and a truncated read sequence. MergeBamAlignment restores this chromosome 10 alignment with a full read sequence and adjusts the CIGAR string to 60M91S, which soft-clips the previously hard-clipped region without loss of alignment specificity. </li>
<li>Both chromosome 2 and chromosome 10 alignments have zero mapping qualities to indicate multiple equally likely mappings. The similar alignment scores of 64 and 50, given by the <code>AS</code> tag, contribute in part to this ambiguity. Additionally, because we asked the aligner to flag shorter split reads as secondary, with the <code>-M</code> option, it assigned a  <code>417</code> <a href="https://broadinstitute.github.io/picard/explain-flags.html">flag</a> to the shorter split chromosome 10 alignment. This makes the chromosome 2 alignment for this read the primary alignment. We set our <code>PRIMARY_ALIGNMENT_STRATEGY</code> to MostDistant which asks the tool to consider the best <em>pair</em> to mark as primary from the primary and secondary records. MergeBamAlignment reassigns the chromosome 10 mapping as the primary alignment (<code>163</code> flag) and the chromosome 2 mapping as secondary (<code>409</code> flag). </li>
<li>MergeBamAlignment updates read group <code>RG</code> information, program group <code>PG</code> information and mate CIGAR <code>MC</code> tags as specified by our command for reads and in the header section. The tool retains <code>SA</code>, <code>MD</code>, <code>NM</code> and <code>AS</code> tags from the alignment, given these are not present in the uBAM. The tool additionally adds <code>UQ</code> (the Phred likelihood of the segment) and <code>MQ</code> (mapping quality of the mate/next segment) tags if applicable. The following table summarizes changes to our tutorial data's tags during the workflow.</li>
</ul>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: center;">original</th>
<th style="text-align: center;">RevertSam</th>
<th style="text-align: center;">BWA MEM</th>
<th style="text-align: center;">MergeBamAlignment</th>
<th></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: center;">RG</td>
<td style="text-align: center;">RG</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">RG</td>
<td>read group</td>
</tr>
<tr>
<td style="text-align: center;">PG</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">PG</td>
<td style="text-align: center;">PG</td>
<td>program group</td>
</tr>
<tr>
<td style="text-align: center;">OC</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td>original cigar</td>
</tr>
<tr>
<td style="text-align: center;">XN</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td># of ambiguous bases in ref</td>
</tr>
<tr>
<td style="text-align: center;">OP</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td>original mapping position</td>
</tr>
<tr>
<td style="text-align: center;">SA</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">SA</td>
<td style="text-align: center;">SA</td>
<td>chimeric alignment</td>
</tr>
<tr>
<td style="text-align: center;">MD</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">MD</td>
<td style="text-align: center;">MD</td>
<td>string for mismatching positions</td>
</tr>
<tr>
<td style="text-align: center;">NM</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">NM</td>
<td style="text-align: center;">NM</td>
<td># of mismatches</td>
</tr>
<tr>
<td style="text-align: center;">AS</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">AS</td>
<td style="text-align: center;">AS</td>
<td>alignment score</td>
</tr>
<tr>
<td style="text-align: center;">UQ</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">UQ</td>
<td>Phred likelihood of the segment</td>
</tr>
<tr>
<td style="text-align: center;">MC</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">MC</td>
<td>CIGAR string for mate</td>
</tr>
<tr>
<td style="text-align: center;">MQ</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">MQ</td>
<td>mapping quality of the mate</td>
</tr>
<tr>
<td style="text-align: center;">OQ</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td>original base quality</td>
</tr>
<tr>
<td style="text-align: center;">XT</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td>tool specific</td>
</tr>
<tr>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">XS</td>
<td style="text-align: center;">XS</td>
<td>BWA's secondary alignment score</td>
</tr>
</tbody>
</table>
<ul>
<li>In our example, we retained the aligner generated <code>XS</code> tag, for secondary alignment scores, with the <code>ATTRIBUTES_TO_RETAIN</code> option.</li>
</ul>
<p>After merging our whole tutorial file, our unmapped read records increases by 620, from 5,334,323 to 5,334,943 due to changes in flag designations and not because any reads failed to map. Our total read records remains the same at 828,846,200 for our 819,728,254 original reads, giving ~1.11% multi-record reads.</p>
<p><a href="#top">back to top</a></p>
<hr />