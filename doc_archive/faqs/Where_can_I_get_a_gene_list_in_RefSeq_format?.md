## Where can I get a gene list in RefSeq format?

http://gatkforums.broadinstitute.org/gatk/discussion/1329/where-can-i-get-a-gene-list-in-refseq-format

<h3>1. About the RefSeq Format</h3>
<p>From the <a href="http://www.ncbi.nlm.nih.gov/refseq/">NCBI RefSeq website</a></p>
<blockquote>
<p>The Reference Sequence (RefSeq) collection aims to provide a comprehensive, integrated, non-redundant, well-annotated set of sequences, including genomic DNA, transcripts, and proteins. RefSeq is a foundation for medical, functional, and diversity studies; they provide a stable reference for genome annotation, gene identification and characterization, mutation and polymorphism analysis (especially RefSeqGene records), expression studies, and comparative analyses.</p>
</blockquote>
<h3>2. In the GATK</h3>
<p>The GATK uses RefSeq in a variety of walkers, from indel calling to variant annotations.  There are many file format flavors of ReqSeq; we've chosen to use the table dump available from the <a href="http://genome.ucsc.edu/cgi-bin/hgTables?command=start">UCSC genome table browser</a>.     </p>
<h3>3. Generating RefSeq files</h3>
<p>Go to the <a href="http://genome.ucsc.edu/cgi-bin/hgTables?command=start">UCSC genome table browser</a>. There are many output options, here are the changes that you'll need to make:   </p>
<pre><code class="pre_md">clade:    Mammal
genome:   Human
assembly: ''choose the appropriate assembly for the reference you're using''
group:    Genes abd Gene Prediction Tracks
track:    RefSeq Genes
table:    refGene
region:   ''choose the genome option''</code class="pre_md"></pre>
<p>Choose a good output filename, something like <code>geneTrack.refSeq</code>, and click the <code>get output</code> button.  You now have your initial RefSeq file, which will not be sorted, and will contain non-standard contigs. To run with the GATK, contigs other than the standard 1-22,X,Y,MT must be removed, and the file sorted in karyotypic order.</p>
<h3>4. Running with the GATK</h3>
<p>You can provide your RefSeq file to the GATK like you would for any other ROD command line argument.  The line would look like the following:</p>
<pre><code class="pre_md">-[arg]:REFSEQ /path/to/refSeq</code class="pre_md"></pre>
<p>Using the filename from above.  </p>
<h4>Warning:</h4>
<p>The GATK automatically adjusts the start and stop position of the records from zero-based half-open intervals (UCSC standard) to one-based closed intervals. </p>
<p>For example: </p>
<pre><code class="pre_md">The first 19 bases in Chromosome one:
Chr1:0-19 (UCSC system)
Chr1:1-19 (GATK)</code class="pre_md"></pre>
<p>All of the GATK output is also in this format, so if you're using other tools or scripts to process RefSeq or GATK output files, you should be aware of this difference.</p>