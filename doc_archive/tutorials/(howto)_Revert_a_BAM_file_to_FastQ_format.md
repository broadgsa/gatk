## (howto) Revert a BAM file to FastQ format

http://gatkforums.broadinstitute.org/gatk/discussion/2908/howto-revert-a-bam-file-to-fastq-format

<h3>NOTE: This tutorial has been replaced by a more recent and much improved version, <a href="http://gatkforums.broadinstitute.org/firecloud/discussion/6484/">Tutorial#6484</a>.</h3>
<hr />
<h4>Objective</h4>
<p>Revert a BAM file back to FastQ. This comes in handy when you receive data that has been processed but not according to GATK Best Practices, and you want to reset and reprocess it properly.</p>
<h4>Prerequisites</h4>
<ul>
<li>Installed HTSlib</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Shuffle the reads in the bam file</li>
<li>Revert the BAM file to FastQ format</li>
<li>Compress the FastQ file </li>
<li>Note for advanced users</li>
</ol>
<hr />
<h3>1. Shuffle the reads in the bam file</h3>
<h4>Action</h4>
<p>Shuffle the reads in the bam file so they are not in a biased order before alignment by running the following HTSlib command: </p>
<pre><code class="pre_md">htscmd bamshuf -uOn 128 aln_reads.bam tmp &gt; shuffled_reads.bam </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a new BAM file containing the original reads, which still retain their mapping information, but now they are no longer sorted. </p>
<p>The aligner uses blocks of paired reads to estimate the insert size. If you don’t shuffle your original bam, the blocks of insert size will not be randomly distributed across the genome, rather they will all come from the same region, biasing the insert size calculation. This is a very important step which is unfortunately often overlooked. </p>
<hr />
<h3>2. Revert the BAM file to FastQ</h3>
<h4>Action</h4>
<p>Revert the BAM file to FastQ format by running the following HTSlib command: </p>
<pre><code class="pre_md">htscmd bam2fq -a shuffled_reads.bam &gt; interleaved_reads.fq </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates an interleaved FastQ file called <code>interleaved_reads.fq</code> containing the now-unmapped paired reads. </p>
<p><em>Interleaved</em> simply means that for each pair of reads in your paired-end data set, both the forward and the reverse reads are in the same file, as opposed to having them in separate files. </p>
<hr />
<h3>3. Compress the FastQ file</h3>
<h4>Action</h4>
<p>Compress the FastQ file to reduce its size using the gzip utility: </p>
<pre><code class="pre_md">gzip interleaved_reads.fq</code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a gzipped FastQ file called <code>interleaved_reads.fq.gz</code>. This file is ready to be used as input for the Best Practices workflow.</p>
<p>BWA handles gzipped fastq files natively, so you don’t need to unzip the file to use it later on. </p>
<hr />
<h3>4. Note for advanced users</h3>
<p>If you’re feeling adventurous, you can do all of the above with this beautiful one-liner, which will save you a heap of time that the program would otherwise spend performing I/O (loading in and writing out data to/from disk): </p>
<pre><code class="pre_md">htscmd bamshuf -uOn 128 aln_reads.bam tmp | htscmd bam2fq -a - | gzip &gt; interleaved_reads.fq.gz </code class="pre_md"></pre>