## (How to) Simulate reads using a reference genome ALT contig

http://gatkforums.broadinstitute.org/gatk/discussion/7859/how-to-simulate-reads-using-a-reference-genome-alt-contig

<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/ca/a90f37bd7d68ff1165411a67c12bd7.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ca/a90f37bd7d68ff1165411a67c12bd7.png" align="right" width="360" style="margin:5px 0px 0px 10px"/></a> This tutorial shows how to generate simulated reads against a specific target sequence. This can be useful, e.g. if you want to simulate reads for an alternate contig in GRCh38/hg38 to see how they end up mapping between the primary assembly versus the alternate contig.</p>
<p>We use external tools to accomplish this. In <strong>Section 1</strong>, we use <a href="http://www.htslib.org/">Samtools</a> to subset the target contig sequence from a reference FASTA file. In <strong>Section 2</strong>, we use <a href="https://github.com/lh3/wgsim">wgsim</a> to generate FASTQ format paired reads against the target contig. The resulting read data is ready for alignment.</p>
<p>This tutorial provides example data for you to follow along and includes a mini-reference FASTA. If you are unfamiliar with terms that describe reference genome components, take a few minutes to study the <em>Dictionary</em> entry <a href="http://gatkforums.broadinstitute.org/dsde/discussion/7857">Reference Genome Components</a>.</p>
<hr />
<h3>Prerequisites and tools involved</h3>
<p>This tutorial uses external tools that may require additional dependencies, e.g. the gcc compiler, that may not be available by default on your system.</p>
<ul>
<li>
<p>After downloading <a href="https://github.com/lh3/wgsim">wgsim</a>, follow instructions to compile. For v0.3.0, the command is as follows. </p>
<pre><code class="pre_md">gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm</code class="pre_md"></pre>
</li>
<li>Samtools. See <a href="http://gatkforums.broadinstitute.org/wdl/discussion/2899">Tutorial#2899</a> for installation instructions.</li>
</ul>
<h3>Download example data</h3>
<ul>
<li>
<p>tutorial_7859.tar.gz (<a href="https://drive.google.com/open?id=0BzI1CyccGsZibnRtQjhaakxobEE">GoogleDrive</a>; <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/tutorials/datasets">ftp site</a>) contains six files. A mini-reference <code>chr19_chr19_KI270866v1_alt.fasta</code> and corresponding <code>.dict</code> dictionary and <code>.fai</code> index, the subset <code>chr19_KI270866v1_alt.fasta</code> and final <code>7859_GPI.read1.fq</code> and <code>7859_GPI.read2.fq</code> FASTQ files.</p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/31/f1f2c77b6efbf9565700516b836914.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/31/f1f2c77b6efbf9565700516b836914.png" align="right" width="360" style="margin:10px 0px 0px 10px"/></a> The mini-reference contains two contigs subset from human GRCh38/hg38: <code>chr19</code> and <code>chr19_KI270866v1_alt</code>. In the tutorial, we simulate reads for the 43,156 bp ALT contig. The ALT contig corresponds to a diverged haplotype of chromosome 19. Specifically, it corresponds to chr19:34350807-34392977, which contains the <em>glucose-6-phosphate isomerase</em> or GPI gene. Part of the ALT contig introduces novel sequence that lacks a corresponding region in the primary assembly.</p>
</li>
</ul>
<hr />
<h2>1. Use Samtools to subset target contig sequence from FASTA reference</h2>
<p>Each contig in the reference FASTA has a header line beginning with <code>&gt;</code> that identifies the contig sequence that follows. We need the exact representation of this header line to subset the target contig sequence. The UNIX command below lists all such headers for the FASTA file.</p>
<pre><code class="pre_md">grep '&gt;' chr19_chr19_KI270866v1_alt.fasta</code class="pre_md"></pre>
<p>This prints the following for our mini-reference <code>chr19_chr19_KI270866v1_alt.fasta</code>.</p>
<pre><code>&gt;chr19
&gt;chr19_KI270866v1_alt</code></pre>
<p>Use the <code>faidx</code> option of Samtools to subset the ALT contig sequence to a new FASTA file, <code>chr19_KI270866v1_alt.fasta</code>.</p>
<pre><code class="pre_md">samtools faidx chr19_chr19_KI270866v1_alt.fasta chr19_KI270866v1_alt &gt; chr19_KI270866v1_alt.fasta</code class="pre_md"></pre>
<hr />
<h3>Optionally introduce variants into reads</h3>
<p>To introduce variants into reads, alter the FASTA sequence at this point before simulating reads. For example, to introduce a simple heterozygous SNP, duplicate the contig information within the file, name the duplicate contig differently, and change the base within the duplicated sequence. Search for the target base's sequence context by using <em>TextEdit</em>'s <em>Find</em> function. Keep in mind FASTA file sequences contain line breaks.  </p>
<p>To generate an alternate FASTA reference based on a VCF of variants, see GATK’s <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_fasta_FastaAlternateReferenceMaker.php">FastaAlternateReferenceMaker</a>.</p>
<hr />
<h2>2. Use wgsim to simulate FASTQ paired reads against the target contig FASTA</h2>
<p>Generate simulated reads from <code>chr19_KI270866v1_alt.fasta</code> with the following command.</p>
<pre><code>wgsim -1151 -2151 -d500 -r0 -e0 -N10000 -R0 -X0 chr19_KI270866v1_alt.fasta 7859_GPI.read1.fq 7859_GPI.read2.fq</code></pre>
<p>This gives two FASTQ files, <code>7859_GPI.read1.fq</code> and <code>7859_GPI.read2.fq</code>, one for each mate of the paired reads.  </p>
<ul>
<li>Each read is 151 bases. Set with <code>-1151</code> and <code>-2151</code> for read1 and read2, respectively.</li>
<li>The outer distance or insert size is 500 bases with a standard deviation of 50. This is set with the <code>-d500</code> parameter.</li>
<li>The files contain 10K read pairs, and this is set by the <code>-N10000</code> parameter.</li>
<li>None of the reads contain indels (<code>-R0</code> &amp; <code>-X0</code>) nor mutations/variants (<code>-r0</code>).</li>
<li>Base quality scales with the value given to <code>-e</code> so we set it to zero (<code>-e0</code>) for base quality scores of <code>I</code>, which is, according to <a href="https://en.wikipedia.org/wiki/FASTQ_format">this page</a> and <a href="http://broadinstitute.github.io/picard/explain-qualities.html">this site</a>, an excellent base quality score equivalent to a Sanger Phred+33 score of 40. </li>
</ul>
<p>For a 43 kb contig, 10K x 2 x 151 reads should give us ~70x hypothetical coverage. Here are two pairs of reads from <code>7859_GPI.read1.fq</code> and <code>7859_GPI.read2.fq</code>.</p>
<p><strong>7859_GPI.read1.fq</strong></p>
<pre><code>@chr19_KI270866v1_alt_40173_40622_0:0:0_0:0:0_0/1
AGGTATGAGGATCTGGGTCTTCCCGTGTCTGAGTAGGTAGCACCTGGCACAGGTATGAGGATATGGGTCTTCCATGTCTGAGGAGGTAGCACCTGGCACAGATATGAGGATCTGCGTCTTCCAGTGTTTGAGGAGGTGAGTTTGGACTCAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@chr19_KI270866v1_alt_30797_31341_0:0:0_0:0:0_1/1
CACCACTGCTGAGCTCAGGCAAGTGCACAAGGAAAGCTGTGGCTCACTGCTCGGCTCCAGCAGAGGTGGTCCCATGGACCACCTGTTGCTACAGAGGGGTCGGCAGCCCTGTCACTCAAGGCAGGGTTTGCTCTGCAAGCTGCCCCAGCCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII</code></pre>
<p><strong>7859_GPI.read2.fq</strong></p>
<pre><code>@chr19_KI270866v1_alt_40173_40622_0:0:0_0:0:0_0/2
AGGGCCAGATCACACCTCCTCAGATATTGACCGACCCAGATCCTTATACCTGCACCAGATCCTACCTCCTCAGGCATTGACAGATCCAGATCCTTATACTTGTGCCAGATCCTACCTCCTTAGACATGGACAGACCCAGATCCTCATACCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@chr19_KI270866v1_alt_30797_31341_0:0:0_0:0:0_1/2
AGGCCCATGAGGTCAGGTCAGTGTTTATTGAGTACCTGCTGCATACCTAGCTTGGGGAAAGGTAGAGAGGCCCTCAGAGAGGCTTGGAGGGCAAGAGCAACCCAGGCAGGATGAGGGCTCCACTTCCACCTGAGGGCGGGCTGAGCTTGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII</code></pre>
<p>All the bases of all the reads from a simulation have the same base quality, and in this instance each base quality is <code>I</code>. Notice the read names of the simulated reads contain useful information, e.g the last read name <code>@chr19_KI270866v1_alt_30797_31341_0:0:0_0:0:0_1/2</code> consists of the following.</p>
<ul>
<li>input FASTA file name <code>chr19_KI270866v1_alt</code></li>
<li>the region that the sequence comes from <code>30797_31341</code></li>
<li>sequencing error, substitutions and gaps <code>0:0:0</code> and the same for the mate <code>0:0:0</code></li>
<li>member pair (0-based indexing in hexadecimal) and mate pair information <code>1/2</code></li>
</ul>
<hr />
<h3>Related resources</h3>
<ul>
<li>To convert FASTQ to BAM, see <strong>Section A</strong> of <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6484">Tutorial#6484</a>.</li>
<li>
<p>To align the reads using BWA-MEM (v0.7.15), you can use the following command. Alternatively, see <a href="https://software.broadinstitute.org/gatk/documentation/article?id=8017">Tutorial#8017</a>.</p>
<pre><code class="pre_md">bwa mem chr19_chr19_KI270866v1_alt.fasta 7859_GPI_alt.read1.fq 7859_GPI_alt.read2.fq &gt; GPI_bwamem.sam</code class="pre_md"></pre>
</li>
<li><a href="https://github.com/nh13/DWGSIM/wiki/Simulating-Reads-with-DWGSIM">DWGSIM Tutorial</a>, a variant of WGSIM</li>
<li><a href="http://gatkforums.broadinstitute.org/dsde/discussion/7857">This Dictionary entry</a> reviews terminology that describes reference genome components.</li>
</ul>
<hr />