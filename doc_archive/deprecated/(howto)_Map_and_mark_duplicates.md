## (howto) Map and mark duplicates

http://gatkforums.broadinstitute.org/gatk/discussion/2799/howto-map-and-mark-duplicates

<blockquote>
<h4>See <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6747">Tutorial#6747</a> for a comparison of <em>MarkDuplicates</em> and <em>MarkDuplicatesWithMateCigar</em>, downloadable example data to follow along, and additional commentary.</h4>
</blockquote>
<hr />
<h4>Objective</h4>
<p>Map the read data to the reference and mark duplicates.</p>
<h4>Prerequisites</h4>
<ul>
<li>This tutorial assumes adapter sequences have been removed.</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Identify read group information</li>
<li>Generate a SAM file containing aligned reads</li>
<li>Convert to BAM, sort and mark duplicates </li>
</ol>
<hr />
<h3>1. Identify read group information</h3>
<p>The read group information is key for downstream GATK functionality. The GATK will not work without a read group tag. Make sure to enter as much metadata as you know about your data in the read group fields provided. For more information about all the possible fields in the @RG tag, take a look at the SAM specification. </p>
<h4>Action</h4>
<p>Compose the read group identifier in the following format:</p>
<pre><code class="pre_md">@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 </code class="pre_md"></pre>
<p>where the <code>\t</code> stands for the tab character. </p>
<hr />
<h3>2. Generate a SAM file containing aligned reads</h3>
<h4>Action</h4>
<p>Run the following BWA command: </p>
<p>In this command, replace read group info by the read group identifier composed in the previous step. </p>
<pre><code class="pre_md">bwa mem -M -R ’&lt;read group info&gt;’ -p reference.fa raw_reads.fq &gt; aligned_reads.sam </code class="pre_md"></pre>
<p>replacing the <code>&lt;read group info&gt;</code> bit with the read group identifier you composed at the previous step. </p>
<p><em>The <code>-M</code> flag causes BWA to mark shorter split hits as secondary (essential for Picard compatibility).</em></p>
<h4>Expected Result</h4>
<p>This creates a file called <code>aligned_reads.sam</code> containing the aligned reads from all input files, combined, annotated and aligned to the same reference. </p>
<p>Note that here we are using a command that is specific for pair end data in an interleaved (read pairs together in the same file, with the forward read followed directly by its paired reverse read) fastq file, which is what we are providing to you as a tutorial file. To map other types of datasets (e.g. single-ended or pair-ended in forward/reverse read files) you will need to adapt the command accordingly. Please see the BWA documentation for exact usage and more options for these commands.</p>
<hr />
<h3>3. Convert to BAM, sort and mark duplicates</h3>
<p>These initial pre-processing operations format the data to suit the requirements of the GATK tools. </p>
<h4>Action</h4>
<p>Run the following Picard command to sort the SAM file and convert it to BAM: </p>
<pre><code class="pre_md">java -jar picard.jar SortSam \ 
    INPUT=aligned_reads.sam \ 
    OUTPUT=sorted_reads.bam \ 
    SORT_ORDER=coordinate </code class="pre_md"></pre>
<h4>Expected Results</h4>
<p>This creates a file called <code>sorted_reads.bam</code> containing the aligned reads sorted by coordinate. </p>
<h4>Action</h4>
<p>Run the following Picard command to mark duplicates: </p>
<pre><code class="pre_md">java -jar picard.jar MarkDuplicates \ 
    INPUT=sorted_reads.bam \ 
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt</code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a sorted BAM file called <code>dedup_reads.bam</code> with the same content as the input file, except that any duplicate reads are marked as such. It also produces a metrics file called <code>metrics.txt</code> containing (can you guess?) metrics.</p>
<h4>Action</h4>
<p>Run the following Picard command to index the BAM file: </p>
<pre><code class="pre_md">java -jar picard.jar BuildBamIndex \ 
    INPUT=dedup_reads.bam </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates an index file for the BAM file called <code>dedup_reads.bai</code>.</p>