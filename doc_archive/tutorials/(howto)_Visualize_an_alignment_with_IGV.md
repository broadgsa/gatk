## (howto) Visualize an alignment with IGV

http://gatkforums.broadinstitute.org/gatk/discussion/6491/howto-visualize-an-alignment-with-igv

<p><a name="top"></a></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c4/d9900bc2df5ef6ed48c9709bb16de9.png" height="270"align="right" border="9"/>
<p>Visualize sequence read alignment data (BAM or SAM) on IGV using this quick-start tutorial. The <em>Integrative Genomics Viewer</em> is a non-GATK tool developed at the Broad Institute that allows for interactive exploration of large genomic datasets. </p>
<h4>Tools involved</h4>
<ul>
<li><a href="https://www.broadinstitute.org/igv/download">IGV downloaded</a> to your desktop </li>
</ul>
<h4>Prerequisites</h4>
<ul>
<li>Coordinate-sorted and aligned BAM or SAM file</li>
<li>Corresponding BAI index </li>
<li>Matching reference genome to which the reads align. See <a href="https://www.broadinstitute.org/igv/Genomes">IGV hosted genomes</a> to check if IGV hosts a reference genome or <a href="https://www.broadinstitute.org/software/igv/LoadGenome">this page</a> for instructions on loading a <code>.genome</code> or <code>FASTA</code> file genome.</li>
</ul>
<h4>Download example data</h4>
<ul>
<li><a href="https://drive.google.com/file/d/0BzI1CyccGsZiYm14LU9YZmhUVzg/view?usp=sharing">tutorial_6491.tar.gz</a> contains a coordinated-sorted BAM and corresponding BAI. Most reads align to a 1 Mbp genomic interval on chromosome 10 (10:96,000,000–97,000,000) of the human GRCh37 reference assembly. Specifically, reads align to GATK bundle's <code>human_g1k_v37_decoy.fasta</code> that corresponds to the <code>Human (1kg, b37+decoy)</code> reference hosted by IGV. </li>
</ul>
<h4>Related resources</h4>
<ul>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/2909/">Tutorial#2909</a> for instructions on coordinate-sorting and indexing alignment data.</li>
<li>See the <a href="https://www.broadinstitute.org/igv/">IGV website</a> for downloads and the extensive <a href="http://www.broadinstitute.org/software/igv/home">user guide</a>. For GATK users, we recommend the sections on <a href="http://www.broadinstitute.org/software/igv/AlignmentData">Viewing Alignments</a> and <a href="http://www.broadinstitute.org/software/igv/viewing_vcf_files">Viewing VCF files</a>.</li>
<li>This <em>How to</em> and its example data are referenced in a larger workflow on <a href="http://gatkforums.broadinstitute.org/discussion/6483">(How to) Efficiently map and clean up short read sequence data</a>. </li>
</ul>
<hr />
<h2>View aligned reads using IGV</h2>
<p>To view aligned reads using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>, the SAM or BAM file must be coordinate-sorted and indexed. </p>
<ol>
<li>Always load the reference genome first. Go to <em>Genomes</em>&gt;<em>Load Genome From Server</em> or load from the drop-down menu in the upper left corner. Select <code>Human (1kg, b37+decoy)</code>.</li>
<li>Load the data file. Go to <em>File</em>&gt;<em>Load from File</em> and select <code>6491_snippet.bam</code>. IGV automatically uses the corresponding <code>6491_snippet.bai</code> index in the same folder.</li>
<li>Zoom in to see alignments. For our tutorial data, copy and paste <code>10:96,867,400-96,869,400</code> into the textbox at the top and press <em>Go</em>. A 2 kbp region of chromosome 10 comes into view as shown in the screenshot above.</li>
</ol>
<p>Alongside read data, IGV automatically generates a coverage track that sums the depth of reads for each genomic position.</p>
<h2>Find a specific read and view as pairs</h2>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/46/7099cda7db1c1281fa9c7078982e78.png" height="270"align="right" border="9" />
<ol>
<li>Right-click on the alignment track and <em>Select by name</em>. Copy and paste <code>H0164ALXX140820:2:2107:7323:30703</code> into the read name textbox and press <em>OK</em>. IGV will highlight two reads corresponding to this query name in bold red. </li>
<li>Right-click on the alignment track and select <em>View as pairs</em>. The two highlighted reads will display in the same row connected by a line as shown in the screenshot.</li>
</ol>
<p>Because IGV holds in memory a limited set of data overlapping with the genomic interval in view (this is what makes IGV fast), the <em>select by name</em> feature also applies only to the data that you call into view. For example, we know this read has a secondary alignment on contig hs37d5 (<code>hs37d5:10,198,000-10,200,000</code>). </p>
<blockquote>
<p>If you jump to this new region, is the read also highlighted in red?</p>
</blockquote>
<hr />
<h3>Some tips</h3>
<p><strong>If you find IGV sluggish</strong>, download a Java Web Start <code>jnlp</code> version of IGV that allows more memory. The highest memory setting as of this writing is 10 GB (RAM) for machines with 64-bit Java. For the tutorial example data, the typical 2 GB allocation is sufficient.</p>
<ul>
<li>To run the <code>jnlp</code> version of IGV, you may need to adjust your system's <em>Java Control Panel</em> settings, e.g. enable Java content in the browser. Also, when first opening the <code>jnlp</code>, overcome Mac OS X's gatekeeper function by right-clicking the saved <code>jnlp</code> and selecting <em>Open</em> <em>with Java Web Start</em>. </li>
</ul>
<p><strong>To change display settings</strong>, check out either the <a href="http://www.broadinstitute.org/software/igv/Preferences#Alignments">Alignment Preferences panel</a> or the <a href="http://www.broadinstitute.org/software/igv/prefs.properties">Alignment track Pop-up menu</a>. For persistent changes to your IGV display settings, use the Preferences panel. For track-by-track changes, use the Pop-up menus.</p>
<p>Default Alignment Preferences settings are tuned to genomic sequence libraries. Go to <em>View</em>&gt;<em>Preferences</em> and make sure the settings under the <em>Alignments</em> tab allows you to view reads of interest, e.g. duplicate reads. </p>
<ul>
<li>IGV saves any changes you make to these settings and applies them to future sessions.</li>
<li>Some changes apply only to new sessions started after the change. </li>
<li>To restore default preferences, delete or rename the <code>prefs.properties</code> file within your system's <code>igv</code> folder. IGV automatically generates a new <code>prefs.properties</code> file with default settings. See [IGV's user guide] for details. </li>
</ul>
<p>After loading data, adjust viewing modes specific to track type by right-clicking on a track to pop up a menu of options. For alignment tracks, these options are described <a href="http://www.broadinstitute.org/software/igv/PopupMenus#AlignmentTrack">here</a>. </p>
<hr />