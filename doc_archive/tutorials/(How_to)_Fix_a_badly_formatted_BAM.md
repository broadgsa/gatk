## (How to) Fix a badly formatted BAM

http://gatkforums.broadinstitute.org/gatk/discussion/2909/how-to-fix-a-badly-formatted-bam

<p><a name="top"></a></p>
<p>Fix a BAM that is not indexed or not sorted, has not had duplicates marked, or is lacking read group information. The options on this page are listed in order of decreasing complexity.</p>
<p>You may ask, is all of this really necessary? The GATK imposes strict formatting guidelines, including requiring certain <a href="http://gatkforums.broadinstitute.org/discussion/6472/">read group information</a>, that other software packages do not require. Although this represents a small additional processing burden upfront, the downstream benefits are numerous, including the ability to process library data individually, and significant gains in speed and parallelization options. </p>
<h4>Prerequisites</h4>
<ul>
<li>Installed Picard tools</li>
<li>If indexing or marking duplicates, then coordinate sorted reads </li>
<li>If coordinate sorting, then reference aligned reads </li>
<li>For each read group ID, a single BAM file. If you have a multiplexed file, separate to individual files per sample. </li>
</ul>
<h4>Jump to a section on this page</h4>
<ol>
<li><a href="#addRG">Add read groups, coordinate sort and index</a> using AddOrReplaceReadGroups</li>
<li><a href="#sort">Coordinate sort and index</a> using SortSam</li>
<li><a href="#index">Index an already coordinate-sorted BAM</a> using BuildBamIndex</li>
<li><a href="#markduplicates">Mark duplicates</a> using MarkDuplicates</li>
</ol>
<h4>Tools involved</h4>
<ul>
<li><a href="http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a></li>
<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#SortSam">SortSam</a></li>
<li><a href="broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex">BuildBamIndex</a></li>
<li><a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates">MarkDuplicates</a></li>
</ul>
<h4>Related resources</h4>
<ul>
<li>Our <a href="http://gatkforums.broadinstitute.org/discussion/6472/">dictionary entry on read groups</a> discusses the importance of assigning appropriate read group fields that differentiate samples and factors that contribute to batch effects, e.g. flow cell lane. Be sure that your read group fields conform to the recommendations.</li>
<li><a href="http://broadinstitute.github.io/picard/command-line-overview.html#Overview">Picard's standard options</a> includes adding <code>CREATE_INDEX</code> to the commands of any of its tools that produce coordinate sorted BAMs.</li>
</ul>
<p><a name="addRG"></a></p>
<hr />
<h3>1. Add read groups, coordinate sort and index using AddOrReplaceReadGroups</h3>
<p>Use Picard's <a href="http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a> to appropriately label read group (<code>@RG</code>) fields, coordinate sort and index a BAM file. Only the five required <code>@RG</code> fields are included in the command shown. Consider the other optional <code>@RG</code> fields for better record keeping. </p>
<pre><code class="pre_md">java -jar picard.jar AddOrReplaceReadGroups \ 
    INPUT=reads.bam \ 
    OUTPUT=reads_addRG.bam \ 
    RGID=H0164.2 \ #be sure to change from default of 1
    RGLB= library1 \ 
    RGPL=illumina \ 
    RGPU=H0164ALXX140820.2 \ 
    RGSM=sample1 \ </code class="pre_md"></pre>
<p>This creates a file called <code>reads_addRG.bam</code> with the same content and sorting as the input file, except the SAM record header's <code>@RG</code> line will be updated with the new information for the specified fields and each read will now have an RG tag filled with the <code>@RG</code> ID field information. Because of this repetition, the length of the <code>@RG</code> ID field contributes to file size.</p>
<p>To additionally coordinate sort by genomic location and create a <code>.bai</code> index, add the following options to the command.</p>
<pre><code class="pre_md">    SORT_ORDER=coordinate \ 
    CREATE_INDEX=true</code class="pre_md"></pre>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee #sets environmental variable for temporary directory</code class="pre_md"></pre>
<p><a name="sort"></a></p>
<hr />
<h3>2. Coordinate sort and index using SortSam</h3>
<p>Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#SortSam">SortSam</a> both sorts and indexes and converts between SAM and BAM formats. For coordinate sorting, reads must be aligned to a reference genome.</p>
<pre><code class="pre_md">java -jar picard.jar SortSam \ 
    INPUT=reads.bam \ 
    OUTPUT=reads_sorted.bam \ 
    SORT_ORDER=coordinate \</code class="pre_md"></pre>
<p>Concurrently index by tacking on the following option.</p>
<pre><code class="pre_md">    CREATE_INDEX=true</code class="pre_md"></pre>
<p>This creates a file called <code>reads_sorted.bam</code> containing reads sorted by genomic location, aka coordinate, and a <code>.bai</code> index file with the same prefix as the output, e.g. <code>reads_sorted.bai</code>, within the same directory.</p>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee #sets environmental variable for temporary directory</code class="pre_md"></pre>
<p><a name="index"></a></p>
<hr />
<h3>3. Index an already coordinate-sorted BAM using BuildBamIndex</h3>
<p>Picard's <a href="broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex">BuildBamIndex</a> allows you to index a BAM that is already coordinate sorted.</p>
<pre><code class="pre_md">java -jar picard.jar BuildBamIndex \ 
    INPUT=reads_sorted.bam </code class="pre_md"></pre>
<p>This creates a <code>.bai</code> index file with the same prefix as the input file, e.g. <code>reads_sorted.bai</code>, within the same directory as the input file. You want to keep this default behavior as many tools require the same prefix and directory location for the pair of files. Note that Picard tools do not systematically create an index file when they output a new BAM file, whereas GATK tools will always output indexed files.</p>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee #sets environmental variable for temporary directory</code class="pre_md"></pre>
<p><a name="markduplicates"></a></p>
<hr />
<h3>4. Mark duplicates using MarkDuplicates</h3>
<p>Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates">MarkDuplicates</a> flags both PCR and optical duplicate reads with a 1024 (0x400) <a href="https://broadinstitute.github.io/picard/explain-flags.html">SAM flag</a>. The input BAM must be coordinate sorted.</p>
<pre><code class="pre_md">java -jar picard.jar MarkDuplicates \ 
    INPUT=reads_sorted.bam \ 
    OUTPUT=reads_markdup.bam \
    METRICS_FILE=metrics.txt \
    CREATE_INDEX=true</code class="pre_md"></pre>
<p>This creates a file called <code>reads_markdup.bam</code> with duplicate reads marked. It also creates a file called <code>metrics.txt</code> containing duplicate read data metrics and a <code>.bai</code> index file.</p>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee #sets environmental variable for temporary directory</code class="pre_md"></pre>
<ul>
<li>During sequencing, which involves PCR amplification within the sequencing machine, by a stochastic process we end up sequencing a proportion of DNA molecules that arose from the same parent insert. To be stringent in our variant discovery, GATK tools discount the duplicate reads as evidence for or against a putative variant. </li>
<li>Marking duplicates is less relevant to targeted amplicon sequencing and RNA-Seq analysis. </li>
<li>Optical duplicates arise from a read being sequenced twice as neighboring clusters.</li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />