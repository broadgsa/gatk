## (How to) Map and clean up short read sequence data efficiently

http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently

<p><a name="top"></a>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/78/7ead8e153ecde1a7e68e8155e0aa7d.png" height="330"align="right" border="30"/>If you are interested in emulating the methods used by the Broad Genomics Platform to pre-process your short read sequencing data, you have landed on the right page. The parsimonious operating procedures outlined in this three-step workflow both maximize data quality, storage and processing efficiency to produce a mapped and <em>clean</em> BAM. This clean BAM is ready for analysis workflows that start with MarkDuplicates. </p>
<p>Since your sequencing data could be in a number of formats, the <strong>first step</strong> of this workflow refers you to specific methods to generate a compatible unmapped BAM (uBAM, <a href="http://gatkforums.broadinstitute.org/discussion/6484/#latest#top">Tutorial#6484</a>) or (uBAM<sup>XT</sup>, <a href="http://gatkforums.broadinstitute.org/discussion/6570/#latest#top">Tutorial#6570 coming soon</a>). Not all unmapped BAMs are equal and these methods emphasize cleaning up prior meta information while giving you the opportunity to assign proper <a href="http://gatkforums.broadinstitute.org/discussion/6472/read-groups">read group fields</a>. The <strong>second step</strong> of the workflow has you marking adapter sequences, e.g. arising from read-through of short inserts, using MarkIlluminaAdapters such that they contribute minimally to alignments and allow the aligner to map otherwise unmappable reads. The <strong>third step</strong> pipes three processes to produce the final BAM. Piping SamToFastq, BWA-MEM and MergeBamAlignment saves time and allows you to bypass storage of larger intermediate FASTQ and SAM files. In particular, MergeBamAlignment merges defined information from the aligned SAM with that of the uBAM to conserve read data, and importantly, it generates additional meta information and unifies meta data. The resulting clean BAM is coordinate sorted, indexed. </p>
<blockquote>
<p>The workflow reflects a <em>lossless</em> operating procedure that retains original sequencing read information within the final BAM file such that data is amenable to reversion and analysis by different means. These practices make scaling up and long-term storage efficient, as one needs only keep the final BAM file.  </p>
</blockquote>
<p><a href="http://gatkforums.broadinstitute.org/profile/Geraldine_VdAuwera">Geraldine_VdAuwera</a> points out that there are many different ways of correctly preprocessing HTS data for variant discovery and ours is only one approach. So keep this in mind.</p>
<p>We present this workflow using real data from a public sample. The original data file, called <code>Solexa-272222</code>, is large at 150 GB. The file contains 151 bp paired PCR-free reads giving 30x coverage of a human whole genome sample referred to as NA12878. The entire sample library was sequenced in a single flow cell lane and thereby assigns all the reads the same read group ID. The example commands work both on this large file and on smaller files containing a subset of the reads, collectively referred to as <code>snippet</code>. NA12878 has a variant in exon 5 of the CYP2C19 gene, on the portion of chromosome 10 covered by the snippet, resulting in a nonfunctional protein. Consistent with GATK's recommendation of using the most up-to-date tools, for the given example results, with the exception of BWA, we used the most current versions of tools as of their testing (September to December 2015). We provide illustrative example results, some of which were derived from processing the original large file and some of which show intermediate stages skipped by this workflow.</p>
<blockquote>
<p>Download example snippet data to follow along the tutorial. </p>
</blockquote>
<p>We welcome feedback. Share your suggestions in the <a href="#bottom">Comments section</a> at the bottom of this page. </p>
<hr />
<h4>Jump to a section</h4>
<ol>
<li><a href="#step1">Generate an unmapped BAM from FASTQ, aligned BAM or BCL</a> </li>
<li><a href="#step2">Mark adapter sequences using MarkIlluminaAdapters</a></li>
<li><a href="#step3">Align reads with BWA-MEM and merge with uBAM using MergeBamAlignment</a>
A. <a href="#step3A">Convert BAM to FASTQ and discount adapter sequences using SamToFastq</a>
B. <a href="#step3B">Align reads and flag secondary hits using BWA-MEM</a>
C. <a href="#step3C">Restore altered data and apply &amp; adjust meta information using MergeBamAlignment</a>
D. <a href="#step3D">Pipe SamToFastq, BWA-MEM and MergeBamAlignment to generate a clean BAM</a></li>
</ol>
<h4>Tools involved</h4>
<ul>
<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters">MarkIlluminaAdapters</a></li>
<li><a href="https://en.wikipedia.org/wiki/Pipeline_(Unix)">Unix pipelines</a></li>
<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq">SamToFastq</a></li>
<li>BWA-MEM (<a href="http://arxiv.org/abs/1303.3997">Li 2013 reference</a>; <a href="http://bioinformatics.oxfordjournals.org/content/30/20/2843.long">Li 2014 benchmarks</a>; <a href="http://bio-bwa.sourceforge.net/">homepage</a>; <a href="http://bio-bwa.sourceforge.net/bwa.shtml">manual</a>)</li>
<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a></li>
</ul>
<h4>Prerequisites</h4>
<ul>
<li>Installed Picard tools</li>
<li>Installed GATK tools</li>
<li>Installed BWA</li>
<li>Reference genome</li>
<li>Illumina or similar tech DNA sequence reads file containing data corresponding to one read group ID. That is, the file contains data from one sample and from one flow cell lane.</li>
</ul>
<h4>Download example data</h4>
<ul>
<li>To download the reference, open <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/">ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/</a> in your browser. Leave the password field blank. Download the following three files (~860 MB) to the same folder: <code>human_g1k_v37_decoy.fasta.gz</code>, <code>.fasta.fai.gz</code>, and <code>.dict.gz</code>. This same reference is available to load in IGV. </li>
<li>I divided the example data into two tarballs: <a href="https://drive.google.com/open?id=0BzI1CyccGsZiTE03V0VDT2VnbFk
">tutorial_6483_piped.tar.gz</a> contains the files for the piped process and <a href="https://drive.google.com/open?id=0BzI1CyccGsZiczdFbTlWSWJjd2M">tutorial_6483_intermediate_files.tar.gz</a> contains the intermediate files produced by running each process independently. The data contain reads originally aligning to a one Mbp genomic interval (10:96,000,000-97,000,000) of GRCh37. The table shows the steps of the workflow, corresponding input and output example data files and approximate minutes and disk space needed to process each step. Additionally, we tabulate the time and minimum storage needed to complete the workflow as presented (piped) or without piping.</li>
</ul>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/b1/7ce5332b0ce9066cac8ac633ddaa00.png" height="150"/>
<h4>Related resources</h4>
<ul>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/2909/">this tutorial</a> to add or replace read groups or coordinate-sort and index a BAM.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/6491/">this tutorial</a> for basic instructions on using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>. </li>
<li>For collecting alignment summary metrics, see <a href="https://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics">CollectAlignmentSummaryMetrics</a>, <a href="http://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics">CollectWgsMetrics</a> and <a href="http://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics">CollectInsertSizeMetrics</a>. See <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html">Picard for metrics definitions</a>. </li>
<li>See <a href="https://broadinstitute.github.io/picard/explain-flags.html">SAM flags</a> to interpret SAM flag values.</li>
<li><a href="http://gatkforums.broadinstitute.org/discussion/2799/#latest">Tutorial#2799</a> gives an example command to mark duplicates.</li>
</ul>
<h4>Other notes</h4>
<ul>
<li>When transforming data files, we stick to using Picard tools over other tools to avoid subtle incompatibilities. </li>
<li>
<p>For large files, (1) use the Java <code>-Xmx</code> setting and (2) set the environmental variable <code>TMP_DIR</code> for a temporary directory.    </p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar MarkIlluminaAdapters \
    TMP_DIR=/path/shlee </code class="pre_md"></pre>
<p>In the command, the <code>-Xmx8G</code> Java option caps the maximum heap size, or memory usage, to eight gigabytes. The path given by <code>TMP_DIR</code> points the tool to scratch space that it can use. These options allow the tool to run without slowing down as well as run without causing an <em>out of memory</em> error. The <code>-Xmx</code> settings we provide here are more than sufficient for most cases. For GATK, 4G is standard, while for Picard less is needed. Some tools, e.g. MarkDuplicates, may require more. These options can be omitted for small files such as the example data and the equivalent command is as follows.</p>
<pre><code class="pre_md">java -jar /path/picard.jar MarkIlluminaAdapters </code class="pre_md"></pre>
<p>To find a system's default maximum heap size, type <code>java -XX:+PrintFlagsFinal -version</code>, and look for <code>MaxHeapSize</code>. Note that any setting beyond available memory spills to storage and slows a system down. If <a href="https://www.broadinstitute.org/gatk/guide/article?id=1975">multithreading</a>, increase memory proportionately to the number of threads. e.g. if 1G is the minimum required for one thread, then use 2G for two threads.</p>
</li>
<li>When I call default options within a command, follow suit to ensure the same results.  </li>
</ul>
<hr />
<p><a name="step1"></a></p>
<h2>1. Generate an unmapped BAM from FASTQ, aligned BAM or BCL</h2>
<p>If you have raw reads data in BAM format with appropriately assigned read group fields, then you can start with step 2. Namely, besides differentiating samples, the read group ID should differentiate factors contributing to technical batch effects, i.e. flow cell lane. If not, you need to reassign read group fields. This <a href="http://gatkforums.broadinstitute.org/discussion/6472/read-groups#latest">dictionary post</a> describes factors to consider and <a href="http://gatkforums.broadinstitute.org/discussion/3060/">this post</a> and <a href="http://gatkforums.broadinstitute.org/discussion/6057/i-have-multiple-read-groups-for-1-sample-how-should-i-pre-process-them#latest">this post</a> provide some strategic advice on handling multiplexed data.  </p>
<ul>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/2909/">this tutorial</a> to add or replace read groups.</li>
</ul>
<p>If your reads are mapped, or in BCL or FASTQ format, then generate an unmapped BAM according to the following instructions.</p>
<ul>
<li>To convert FASTQ or revert aligned BAM files, follow directions in <a href="http://gatkforums.broadinstitute.org/discussion/6484/#latest#top">Tutorial#6484</a>. The resulting uBAM needs to have its adapter sequences marked as outlined in the next step (step 2).</li>
<li>To convert an Illumina Base Call files (BCL) use <a href="http://broadinstitute.github.io/picard/command-line-overview.html#IlluminaBasecallsToSam">IlluminaBasecallsToSam</a>. The tool marks adapter sequences at the same time. The resulting uBAM<sup>XT</sup> has adapter sequences marked with the XT tag so you can skip step 2 of this workflow and go directly to step 3. The corresponding <a href="http://gatkforums.broadinstitute.org/discussion/6570/">Tutorial#6570</a> is coming soon.</li>
</ul>
<blockquote>
<p>See if you can revert <code>6483_snippet.bam</code>, containing 279,534 aligned reads, to the unmapped <code>6383_snippet_revertsam.bam</code>, containing 275,546 reads. </p>
</blockquote>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step2"></a></p>
<h2>2. Mark adapter sequences using MarkIlluminaAdapters</h2>
<p><a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkIlluminaAdapters">MarkIlluminaAdapters</a> adds the XT tag to a read record to mark the 5' start position of the specified adapter sequence and produces a metrics file. Some of the marked adapters come from concatenated adapters that randomly arise from the primordial soup that is a PCR reaction. Others represent read-through to 3' adapter ends of reads and arise from insert sizes that are shorter than the read length. In some instances read-though can affect the majority of reads in a sample, e.g. in Nextera library samples over-titrated with transposomes, and render these reads unmappable by certain aligners. Tools such as SamToFastq use the XT tag in various ways to effectively remove adapter sequence contribution to read alignment and alignment scoring metrics. Depending on your library preparation, insert size distribution and read length, expect varying amounts of such marked reads. </p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar MarkIlluminaAdapters \
I=6483_snippet_revertsam.bam \
O=6483_snippet_markilluminaadapters.bam \
M=6483_snippet_markilluminaadapters_metrics.txt \ #naming required
TMP_DIR=/path/shlee #optional to process large files</code class="pre_md"></pre>
<p>This produces two files. (1) The metrics file, <code>6483_snippet_markilluminaadapters_metrics.txt</code> bins the number of tagged adapter bases versus the number of reads. (2) The <code>6483_snippet_markilluminaadapters.bam</code> file is identical to the input BAM, <code>6483_snippet_revertsam.bam</code>, except reads with adapter sequences will be marked with a tag in XT:i:# format, where # denotes the 5' starting position of the adapter sequence. At least six bases are required to mark a sequence. Reads without adapter sequence remain untagged.  </p>
<ul>
<li>By default, the tool uses Illumina adapter sequences. This is sufficient for our example data. </li>
<li>Adjust the default standard Illumina adapter sequences to any adapter sequence using the <code>FIVE_PRIME_ADAPTER</code> and <code>THREE_PRIME_ADAPTER</code> parameters. To clear and add new adapter sequences first set <code>ADAPTERS</code> to 'null' then specify each sequence with the parameter.  </li>
</ul>
<p>We plot the metrics data that is in <a href="https://www.broadinstitute.org/gatk/guide/article?id=1244">GATKReport file format</a> using RStudio, and as you can see, marked bases vary in size up to the full length of reads.
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/fd/a8015ccd957ab6be075ff3609c8896.png" height="230" border="9" /> <img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a1/59e1837963fe37d0577d466f3c56b2.png"height="230" border="9" /></p>
<blockquote>
<p>Do you get the same number of marked reads? <code>6483_snippet</code> marks 448 reads (0.16%) with XT, while the original <code>Solexa-272222</code> marks 3,236,552 reads (0.39%). </p>
</blockquote>
<p>Below, we show a read pair marked with the XT tag by MarkIlluminaAdapters. The insert region sequences for the reads overlap by a length corresponding approximately to the XT tag value. For XT:i:20, the majority of the read is adapter sequence. The same read pair is shown after SamToFastq transformation, where adapter sequence base quality scores have been set to 2 (# symbol), and after MergeBamAlignment, which restores original base quality scores. </p>
<p><strong>Unmapped uBAM (step 1)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/92/efe96619d7c3f637ec07c7844540c3.png" />
<p><strong>After MarkIlluminaAdapters (step 2)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/05/0ddc4c9f0900c65b6d4dbdb3078c28.png" />
<p><strong>After SamToFastq (step 3)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/8e/9b0a614feda8111a5d2cb81badb05d.png" />
<p><strong>After MergeBamAlignment (step 3)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/27/47b6523c1fa5936ed0810709985ed7.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3"></a></p>
<h2>3. Align reads with BWA-MEM and merge with uBAM using MergeBamAlignment</h2>
<p>This step actually pipes three processes, performed by three different tools. Our tutorial example files are small enough to easily view, manipulate and store, so any difference in piped or independent processing will be negligible. For larger data, however, using <a href="https://en.wikipedia.org/wiki/Pipeline_(Unix)">Unix pipelines</a> can add up to significant savings in processing time and storage.  </p>
<blockquote>
<p>Not all tools are amenable to piping and piping the wrong tools or wrong format can result in anomalous data.</p>
</blockquote>
<p>The three tools we pipe are SamToFastq, BWA-MEM and MergeBamAlignment. By piping these we bypass storage of larger intermediate FASTQ and SAM files. We additionally save time by eliminating the need for the processor to read in and write out data for two of the processes, as piping retains data in the processor's input-output (I/O) device for the next process. </p>
<p>To make the information more digestible, we will first talk about each tool separately. At the end of the section, we provide the piped command.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3A"></a></p>
<h3>3A. Convert BAM to FASTQ and discount adapter sequences using SamToFastq</h3>
<p>Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq">SamToFastq</a> takes read identifiers, read sequences, and base quality scores to write a Sanger FASTQ format file. We use additional options to effectively remove previously marked adapter sequences, in this example marked with an XT tag. By specifying <code>CLIPPING_ATTRIBUTE</code>=XT and <code>CLIPPING_ACTION</code>=2, SamToFastq changes the quality scores of bases marked by XT to two--a rather low score in the Phred scale. This effectively removes the adapter portion of sequences from contributing to downstream read alignment and alignment scoring metrics. </p>
<p><strong>Illustration of an intermediate step unused in workflow. See <a href="#step3D">piped command</a>.</strong></p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar SamToFastq \
I=6483_snippet_markilluminaadapters.bam \
FASTQ=6483_snippet_samtofastq_interleaved.fq \
CLIPPING_ATTRIBUTE=XT \
CLIPPING_ACTION=2 \
INTERLEAVE=true \ 
NON_PF=true \
TMP_DIR=/path/shlee #optional to process large files         </code class="pre_md"></pre>
<p>This produces a FASTQ file in which all extant meta data, i.e. read group information, alignment information, flags and tags are purged. What remains are the read query names prefaced with the <code>@</code> symbol, read sequences and read base quality scores. </p>
<ul>
<li>For our paired reads example file we set SamToFastq's <code>INTERLEAVE</code> to true. During the conversion to FASTQ format, the query name of the reads in a pair are marked with /1 or /2 and paired reads are retained in the same FASTQ file. <a href="http://bio-bwa.sourceforge.net/bwa.shtml">BWA aligner</a> accepts interleaved FASTQ files given the <code>-p</code> option. </li>
<li>We change the <code>NON_PF</code>, aka <code>INCLUDE_NON_PF_READS</code>, option from default to true. SamToFastq will then retain reads marked by what <a href="https://github.com/samtools/hts-specs/issues/85">some consider an archaic 0x200 flag bit</a> that denotes reads that do not pass quality controls, aka reads failing platform or vendor quality checks. Our tutorial data do not contain such reads and we call out this option for illustration only.</li>
<li>Other CLIPPING_ACTION options include (1) X to hard-clip, (2) N to change bases to Ns or (3) another number to change the base qualities of those positions to the given value.</li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3B"></a></p>
<h3>3B. Align reads and flag secondary hits using BWA-MEM</h3>
<p>In this workflow, alignment is the most compute intensive and will take the longest time. GATK's variant discovery workflow recommends Burrows-Wheeler Aligner's maximal exact matches (BWA-MEM) algorithm (<a href="http://arxiv.org/abs/1303.3997">Li 2013 reference</a>; <a href="http://bioinformatics.oxfordjournals.org/content/30/20/2843.long">Li 2014 benchmarks</a>; <a href="http://bio-bwa.sourceforge.net/">homepage</a>; <a href="http://bio-bwa.sourceforge.net/bwa.shtml">manual</a>). BWA-MEM is suitable for aligning high-quality long reads ranging from 70 bp to 1 Mbp against a large reference genome such as the human genome.  </p>
<ul>
<li>Aligning our <code>snippet</code> reads against either a portion or the whole genome is not equivalent to aligning our original <code>Solexa-272222</code> file, merging and taking a new <code>slice</code> from the same genomic interval. </li>
<li>For the tutorial, we use BWA v 0.7.7.r441, the same aligner used by the Broad Genomics Platform as of this writing (9/2015).</li>
<li>As mentioned, alignment is a compute intensive process. For faster processing, use a reference genome with decoy sequences, also called a <a href="http://www.cureffi.org/2013/02/01/the-decoy-genome/">decoy genome</a>. For example, the Broad's Genomics Platform uses an Hg19/GRCh37 reference sequence that includes Ebstein-Barr virus (EBV) sequence to soak up reads that fail to align to the human reference that the aligner would otherwise spend an inordinate amount of time trying to align as split reads. <a href="https://www.broadinstitute.org/gatk/guide/article.php?id=1213">GATK's resource bundle</a> provides a standard decoy genome from the <a href="http://www.1000genomes.org/">1000 Genomes Project</a>.</li>
<li>
<p>BWA alignment requires an indexed reference genome file. Indexing is specific to algorithms. To index the human genome for BWA, we apply BWA's <code>index</code> function on the reference genome file, e.g. <code>human_g1k_v37_decoy.fasta</code>. This produces five index files with the extensions <code>amb</code>, <code>ann</code>, <code>bwt</code>, <code>pac</code> and <code>sa</code>. </p>
<pre><code class="pre_md">bwa index -a bwtsw human_g1k_v37_decoy.fasta</code class="pre_md"></pre>
</li>
</ul>
<p>The example command below aligns our example data against the GRCh37 genome. The tool automatically locates the index files within the same folder as the reference FASTA file. </p>
<p><strong>Illustration of an intermediate step unused in workflow. See <a href="#step3D">piped command</a>.</strong></p>
<pre><code class="pre_md">/path/bwa mem -M -t 7 -p /path/human_g1k_v37_decoy.fasta \ 
6483_snippet_samtofastq_interleaved.fq &gt; 6483_snippet_bwa_mem.sam</code class="pre_md"></pre>
<p>This command takes the FASTQ file, <code>6483_snippet_samtofastq_interleaved.fq</code>, and produces an aligned SAM format file, <code>6483_snippet_unthreaded_bwa_mem.sam</code>, containing read alignment information, an automatically generated program group record and reads sorted in the same order as the input FASTQ file. Aligner-assigned alignment information, flag and tag values reflect each read's or split read segment's best sequence match and does not take into consideration whether pairs are mapped optimally or if a mate is unmapped. Added tags include the aligner-specific <code>XS</code> tag that marks secondary alignment scores in XS:i:# format. This tag is given for each read even when the score is zero and even for unmapped reads. The program group record (@PG) in the header gives the program group ID, group name, group version and recapitulates the given command. Reads are sorted by query name. For the given version of BWA, the aligned file is in SAM format even if given a BAM extension. </p>
<blockquote>
<p>Does the aligned file contain read group information? </p>
</blockquote>
<p>We invoke three options in the command. </p>
<ul>
<li><code>-M</code> to flag shorter split hits as secondary.
This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments.  </li>
<li><code>-p</code> to indicate the given file contains interleaved paired reads.</li>
<li>
<p><code>-t</code> followed by a number for the number of processor threads to use concurrently. Here we use seven threads which is one less than the total threads available on my Mac laptap. Check your server or system's total number of threads with the following command provided by <a href="http://gatkforums.broadinstitute.org/profile/KateN">KateN</a>.</p>
<pre><code class="pre_md">getconf _NPROCESSORS_ONLN </code class="pre_md"></pre>
</li>
</ul>
<p>In the example data, all of the 1211 <em>unmapped</em> reads each have an asterisk (*) in column 6 of the SAM record, where a read typically records its CIGAR string. The asterisk represents that the CIGAR string is unavailable. The several asterisked reads I examined are recorded as mapping exactly to the same location as their <em>mapped</em> mates but with MAPQ of zero. Additionally, the asterisked reads had varying noticeable amounts of low base qualities, e.g. strings of #s, that corresponded to original base quality calls and not those changed by SamToFastq. This accounting by BWA allows these pairs to always list together, even when the reads are coordinate-sorted, and leaves a pointer to the genomic mapping of the mate of the unmapped read. For the example read pair shown below, comparing sequences shows no apparent overlap, with the highest identity at 72% over 25 nts.</p>
<p><strong>After MarkIlluminaAdapters (step 2)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/3c/b4f8b576ae39bcfe8f4b2ebddacd78.png" />
<p><strong>After BWA-MEM (step 3)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a6/d9c889d4dbc5dd3c5d8c79a28086d7.png" />
<p><strong>After MergeBamAlignment (step 3)</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/94/2ff9e8345ae113e4a806e1b4745980.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3C"></a></p>
<h3>3C. Restore altered data and apply &amp; adjust meta information using MergeBamAlignment</h3>
<p><a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a> is a beast of a tool, so its introduction is longer. It does more than is implied by its name. Explaining these features requires I fill you in on some background. </p>
<p>Broadly, the tool merges defined information from the unmapped BAM (uBAM, step 1) with that of the aligned BAM (step 3) to conserve read data, e.g. original read information and base quality scores. The tool also generates additional meta information based on the information generated by the aligner, which may alter aligner-generated designations, e.g. mate information and secondary alignment flags. The tool then makes adjustments so that all meta information is congruent, e.g. read and mate strand information based on proper mate designations. We ascribe the resulting BAM as <em>clean</em>. </p>
<p>Specifically, the aligned BAM generated in step 3 lacks read group information and certain tags--the UQ (Phred likelihood of the segment), MC (CIGAR string for mate) and MQ (mapping quality of mate) tags. It has hard-clipped sequences from split reads and altered base qualities. The reads also have what some call mapping artifacts but what are really just features we should not expect from our aligner. For example, the meta information so far does not consider whether pairs are optimally mapped and whether a mate is unmapped (in reality or for accounting purposes). Depending on these assignments, MergeBamAlignment adjusts the read and read mate strand orientations for reads in a proper pair. Finally, the alignment records are sorted by query name. We would like to fix all of these issues before taking our data to a variant discovery workflow.</p>
<p>Enter MergeBamAlignment. As the tool name implies, MergeBamAlignment applies read group information from the uBAM and retains the program group information from the aligned BAM. In restoring original sequences, the tool adjusts CIGAR strings from hard-clipped to soft-clipped. If the alignment file is missing reads present in the unaligned file, then these are retained as unmapped records. Additionally, MergeBamAlignment evaluates primary alignment designations according to a user-specified strategy, e.g. for optimal <em>mate pair</em> mapping, and changes <em>secondary alignment</em> and <em>mate unmapped</em> <a href="https://broadinstitute.github.io/picard/explain-flags.html">flags</a> based on its calculations. Additional for desired congruency. I will soon explain these and additional changes in more detail and show a read record to illustrate.</p>
<blockquote>
<p>Consider what <code>PRIMARY_ALIGNMENT_STRATEGY</code> option best suits your samples. MergeBamAlignment applies this strategy to a read for which the aligner has provided more than one primary alignment, and for which one is designated primary by virtue of another record being marked secondary. MergeBamAlignment considers and switches only existing primary and secondary designations. Therefore, it is critical that these were previously flagged.</p>
</blockquote>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/49/6079b9ddfca33aedd25ae90ef133f4.png" height="240"align="right" border="9"/> A read with multiple alignment records may map to multiple loci or may be chimeric--that is, splits the alignment. It is possible for an aligner to produce multiple alignments as well as multiple primary alignments, e.g. in the case of a linear alignment set of split reads. When one alignment, or alignment set in the case of chimeric read records, is designated primary, others are designated either secondary or supplementary. Invoking the <code>-M</code> option, we had BWA mark the record with the longest aligning section of split reads as primary and all other records as secondary. MergeBamAlignment further adjusts this secondary designation and adds the read mapped in proper pair (0x2) and mate unmapped (0x8) flags. The tool then adjusts the strand orientation flag for a read (0x10) and it proper mate (0x20). </p>
<p>In the command, we change <code>CLIP_ADAPTERS</code>, <code>MAX_INSERTIONS_OR_DELETIONS</code> and <code>PRIMARY_ALIGNMENT_STRATEGY</code> values from default, and invoke other optional parameters. The path to the reference FASTA given by <code>R</code> should also contain the <a href="http://gatkforums.broadinstitute.org/discussion/1601/">corresponding <code>.dict</code> sequence dictionary</a> with the same prefix as the reference FASTA. It is imperative that both the uBAM and aligned BAM are both sorted by queryname.</p>
<p><strong>Illustration of an intermediate step unused in workflow. See <a href="#step3D">piped command</a>.</strong></p>
<pre><code class="pre_md">java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
R=/path/Homo_sapiens_assembly19.fasta \ 
UNMAPPED_BAM=6383_snippet_revertsam.bam \ 
ALIGNED_BAM=6483_snippet_bwa_mem.sam \ #accepts either SAM or BAM
O=6483_snippet_mergebamalignment.bam \
CREATE_INDEX=true \ #standard Picard option for coordinate-sorted outputs
ADD_MATE_CIGAR=true \ #default; adds MC tag
CLIP_ADAPTERS=false \ #changed from default
CLIP_OVERLAPPING_READS=true \ #default; soft-clips ends so mates do not extend past each other
INCLUDE_SECONDARY_ALIGNMENTS=true \ #default
MAX_INSERTIONS_OR_DELETIONS=-1 \ #changed to allow any number of insertions or deletions
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \ #changed from default BestMapq
ATTRIBUTES_TO_RETAIN=XS \ #specify multiple times to retain tags starting with X, Y, or Z 
TMP_DIR=/path/shlee #optional to process large files</code class="pre_md"></pre>
<p>This generates a coordinate-sorted and <em>clean</em> BAM, <code>6483_snippet_mergebamalignment.bam</code>, and corresponding <code>.bai</code> index. These are ready for analyses starting with MarkDuplicates. The two bullet-point lists below describe changes to the resulting file. The first list gives general comments on select parameters and the second describes some of the notable changes to our example data.</p>
<p><strong>Comments on select parameters</strong></p>
<ul>
<li>Setting <code>PRIMARY_ALIGNMENT_STRATEGY</code>to MostDistant marks primary alignments based on the alignment <em>pair</em> with the largest insert size. This strategy is based on the premise that of chimeric sections of a read aligning to consecutive regions, the alignment giving the largest insert size with the mate gives the most information.</li>
<li>It may well be that alignments marked as secondary represent interesting biology, so we retain them with the <code>INCLUDE_SECONDARY_ALIGNMENTS</code> parameter. </li>
<li>Setting <code>MAX_INSERTIONS_OR_DELETIONS</code> to -1 retains reads irregardless of the number of insertions and deletions. The default is 1.</li>
<li>Because we leave the <code>ALIGNER_PROPER_PAIR_FLAGS</code> parameter at the default false value, MergeBamAlignment will reassess and reassign <em>proper pair</em> designations made by the aligner. These are explained below using the example data.</li>
<li><code>ATTRIBUTES_TO_RETAIN</code> is specified to carryover the XS tag from the alignment, which reports BWA-MEM's suboptimal alignment scores. My impression is that this is the next highest score for any alternative or additional alignments BWA considered, whether or not these additional alignments made it into the final aligned records. (<a href="http://www.broadinstitute.org/software/igv/BLAT">IGV's BLAT feature</a> allows you to search for additional sequence matches). For our tutorial data, this is the only additional unaccounted tag from the alignment. The XS tag in unnecessary for the Best Practices Workflow and is not retained by the Broad Genomics Platform's pipeline. We retain it here not only to illustrate that the tool carries over select alignment information only if asked, but also because I think it prudent. Given how compute intensive the alignment process is, the additional ~1% gain in the <code>snippet</code> file size seems a small price against having to rerun the alignment because we realize later that we want the tag. </li>
<li>Setting <code>CLIP_ADAPTERS</code> to false leaves reads unclipped.</li>
<li>By default the merged file is coordinate sorted. We set <code>CREATE_INDEX</code> to true to additionally create the <code>bai</code> index.</li>
<li>We need not invoke <code>PROGRAM</code> options as BWA's program group information is sufficient and is retained in the merging. </li>
<li>As a standalone tool, we would normally feed in a BAM file for <code>ALIGNED_BAM</code> instead of the much larger SAM. We will be piping this step however and so need not add an extra conversion to BAM. </li>
</ul>
<p><strong>Description of changes to our example data</strong></p>
<ul>
<li>MergeBamAlignment merges header information from the two sources that define read groups (@RG) and program groups (@PG) as well as reference contigs. </li>
<li><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/52/f859f64f062395ea60ca3acbda0ff0.png" height="180"align="right" border="9"/>Tags are updated for our example data as shown in the table. The tool retains SA, MD, NM and AS tags from the alignment, given these are not present in the uBAM. The tool additionally adds UQ (the Phred likelihood of the segment), MC (mate CIGAR string) and MQ (mapping quality of the mate/next segment) tags if applicable. For unmapped reads (marked with an <code>*</code> asterisk in column 6 of the SAM record), the tool removes AS and XS tags and assigns MC (if applicable), PG and RG tags. This is illustrated for example read <code>H0164ALXX140820:2:1101:29704:6495</code> in the BWA-MEM section of this document.</li>
<li>Original base quality score restoration is illustrated in <a href="#step2">step 2</a>. </li>
</ul>
<p>The example below shows a read pair for which MergeBamAlignment adjusts multiple information fields, and these changes are described in the remaining bullet points.</p>
<ul>
<li>MergeBamAlignment changes hard-clipping to soft-clipping, e.g. 96H55M to 96S55M, and restores corresponding truncated sequences with the original full-length read sequence. </li>
<li>The tool reorders the read records to reflect the chromosome and contig ordering in the header and the genomic coordinates for each.</li>
<li>MergeBamAlignment's MostDistant <code>PRIMARY_ALIGNMENT_STRATEGY</code> asks the tool to consider the best <em>pair</em> to mark as primary from the primary and secondary records. In this pair, one of the reads has two alignment loci, on <a href="https://wiki.dnanexus.com/Scientific-Notes/human-genome#The-%22b37+decoy%22-/-%22hs37d5%22-extensions-%28by-the-1000-Genomes-Project-Phase-II%29.">contig hs37d5</a> and on chromosome 10. The two loci align 115 and 55 nucleotides, respectively, and the aligned sequences are identical by 55 bases. <a href="https://broadinstitute.github.io/picard/explain-flags.html">Flag values</a> set by BWA-MEM indicate the contig hs37d5 record is primary and the <em>shorter</em> chromosome 10 record is secondary. For this chimeric read, MergeBamAlignment reassigns the chromosome 10 mapping as the primary alignment and the contig hs37d5 mapping as secondary (0x100 flag bit). </li>
<li>In addition, MergeBamAlignment designates each record on chromosome 10 as <em>read mapped in proper pair</em> (0x2 flag bit) and the contig hs37d5 mapping as <em>mate unmapped</em> (0x8 flag bit). <a href="http://www.broadinstitute.org/igv/">IGV</a>'s <em>paired reads mode</em> displays the two chromosome 10 mappings as a pair after these MergeBamAlignment adjustments. </li>
<li>MergeBamAlignment adjusts <em>read reverse strand</em> (0x10 flag bit) and <em>mate reverse strand</em> (0x20 flag bit) flags consistent with changes to the <em>proper pair</em> designation. For our non-stranded DNA-Seq library alignments displayed in IGV, a read pointing rightward is in the forward direction (absence of 0x10 flag) and a read pointing leftward is in the reverse direction (flagged with 0x10). In a typical pair, where the rightward pointing read is to the left of the leftward pointing read, the left read will also have the <em>mate reverse strand</em> (0x20) flag. </li>
</ul>
<blockquote>
<p>Two distinct classes of <em>mate unmapped</em> read records are now present in our example file: (1) reads whose mates truly failed to map and are marked by an asterisk <code>*</code> in column 6 of the SAM record and (2) multimapping reads whose mates are in fact mapped but in a proper pair that excludes the particular read record. Each of these two classes of <em>mate unmapped</em> reads can contain multimapping reads that map to two or more locations. </p>
</blockquote>
<p>Comparing <code>6483_snippet_bwa_mem.sam</code> and <code>6483_snippet_mergebamalignment.bam</code>, we see the number<em>unmapped reads</em> remains the same at 1211, while the number of records with the <em>mate unmapped</em> flag increases by 1359, from 1276 to 2635. These now account for 0.951% of the 276,970 read records.  </p>
<blockquote>
<p>For <code>6483_snippet_mergebamalignment.bam</code>, how many additional unique reads become <em>mate unmapped</em>? </p>
</blockquote>
<p><strong>After BWA-MEM alignment</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/cf/d2bb7a2995b2c829e33ff9540c0d3d.png" />
<p><strong>After MergeBamAlignment</strong></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ba/ffd65bd327898381f8df902efb31fe.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="step3D"></a></p>
<h3>3D. Pipe SamToFastq, BWA-MEM and MergeBamAlignment to generate a clean BAM</h3>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ea/4f93dc6b258935c49b1aa0f8a27a01.jpg" height="120"align="left" border="27"/> We pipe the three tools described above to generate an aligned BAM file sorted by query name. In the piped command, the commands for the three processes are given together, separated by a <a href="https://en.wikipedia.org/wiki/Vertical_bar">vertical bar</a> <code>|</code>. We also replace each intermediate output and input file name with a symbolic path to the system's output and input devices, here <code>/dev/stdout</code> and <code>/dev/stdin</code>, respectively. We need only provide the first input file and name the last output file.</p>
<p>Before using a piped command, we should <a href="https://sipb.mit.edu/doc/safe-shell/">ask UNIX to stop the piped command</a> if any step of the pipe should error and also return to us the error messages. Type the following into your shell to set these UNIX options.</p>
<pre><code class="pre_md">set -o pipefail</code class="pre_md"></pre>
<p><strong>Overview of command structure</strong></p>
<pre><code class="pre_md">[SamToFastq] | [BWA-MEM] | [MergeBamAlignment]</code class="pre_md"></pre>
<p><strong>Piped command</strong></p>
<pre><code class="pre_md">java -Xmx8G -jar /path/picard.jar SamToFastq \
I=6483_snippet_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=/path/shlee | \ 
/path/bwa mem -M -t 7 -p /path/Homo_sapiens_assembly19.fasta /dev/stdin | \  
java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=6383_snippet_revertsam.bam \ 
OUTPUT=6483_snippet_piped.bam \
R=/path/Homo_sapiens_assembly19.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/path/shlee</code class="pre_md"></pre>
<p>The piped output file, <code>6483_snippet_piped.bam</code>, is for all intensive purposes the same as <code>6483_snippet_mergebamalignment.bam</code>, produced by running MergeBamAlignment separately without piping. However, the resulting files, as well as new runs of the workflow on the same data, have the potential to differ in small ways because each uses a different alignment instance. </p>
<blockquote>
<p>How do these small differences arise? </p>
</blockquote>
<p>Counting the number of <em>mate unmapped</em> reads shows that this number remains unchanged for the two described workflows. Two counts emitted at the end of the process updates, that also remain constant for these instances, are the number of alignment records and the number of unmapped reads. </p>
<pre><code class="pre_md">INFO    2015-12-08 17:25:59 AbstractAlignmentMerger Wrote 275759 alignment records and 1211 unmapped reads.</code class="pre_md"></pre>
<p><a href="#top">back to top</a></p>
<hr />
<h3>Some final remarks</h3>
<p>We have produced a <em>clean</em> BAM that is coordinate-sorted and indexed, in an efficient manner that minimizes processing time and storage needs. The file is ready for marking duplicates as outlined in <a href="http://gatkforums.broadinstitute.org/discussion/2799/#latest">Tutorial#2799</a>. Additionally, we can now free up storage on our file system by deleting the original file we started with, the uBAM and the uBAM<sup>XT</sup>. We sleep well at night knowing that the clean BAM retains all original information.</p>
<p>We have two final comments (1) on multiplexed samples and (2) on fitting this workflow into a larger workflow.</p>
<p>For multiplexed samples, first perform the workflow steps on a file representing one sample and one lane. Then mark duplicates. Later, after some steps in the GATK's variant discovery workflow, and after aggregating files from the same sample from across lanes into a single file, mark duplicates again. These two marking steps ensure you find both optical and PCR duplicates.</p>
<p>For workflows that nestle this pipeline, consider additionally optimizing java jar's parameters for SamToFastq and MergeBamAlignment. For example, the following are the additional settings used by the Broad Genomics Platform in the piped command for very large data sets.  </p>
<pre><code class="pre_md">    java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m -jar /path/picard.jar SamToFastq ...

    java -Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m -jar /path/picard.jar MergeBamAlignment ...</code class="pre_md"></pre>
<p>I give my sincere thanks to Julian Hess, the GATK team and the Data Sciences and Data Engineering (DSDE) team members for all their help in writing this and related documents.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="bottom"></a></p>