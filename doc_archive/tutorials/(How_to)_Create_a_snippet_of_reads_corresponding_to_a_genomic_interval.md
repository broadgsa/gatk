## (How to) Create a snippet of reads corresponding to a genomic interval

http://gatkforums.broadinstitute.org/gatk/discussion/6517/how-to-create-a-snippet-of-reads-corresponding-to-a-genomic-interval

<h4>Tools involved</h4>
<ul>
<li><a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php">PrintReads</a></li>
</ul>
<h4>Prerequisites</h4>
<ul>
<li>Installed GATK tools</li>
<li>Reference genome</li>
<li>Coordinate-sorted, aligned and indexed BAM </li>
</ul>
<h4>Download example data</h4>
<ul>
<li>Use the <a href="http://gatkforums.broadinstitute.org/discussion/4610/">advanced tutorial bundle</a>'s human_g1k_v37_decoy.fasta as reference </li>
<li><a href="https://drive.google.com/open?id=0BzI1CyccGsZiTmlDLW13MXdTSG8">tutorial_6517.tar.gz</a> contains four files: 6517_2Mbp_input.bam and .bai covering reads aligning to 10:90,000,000-92,000,000 and 6517_1Mbp_output.bam and .bai covering 10:91,000,000-92,000,000</li>
</ul>
<h4>Related resources</h4>
<ul>
<li>This <em>How to</em> is referenced in a tutorial on <a href="http://gatkforums.broadinstitute.org/discussion/6484/">(How to) Generate an unmapped BAM (uBAM)</a>. </li>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/2909/">this tutorial</a> to coordinate-sort and index a BAM.</li>
<li>See <a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--unsafe">here</a> for command line parameters accepted by all GATK tools.</li>
<li>For applying interval lists, e.g. to whole exome data, see <a href="http://gatkforums.broadinstitute.org/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals">When should I use L to pass in a list of intervals</a>.</li>
</ul>
<hr />
<h3>Create a snippet of reads corresponding to a genomic interval using PrintReads</h3>
<p>PrintReads merges or subsets sequence data. The tool automatically applies MalformedReadFilter and BadCigarFilter to filter out certain types of reads that cause problems for downstream GATK tools, e.g. reads with mismatching numbers of bases and base qualities or reads with CIGAR strings containing the N operator.  </p>
<ul>
<li>To create a test snippet of RNA-Seq data that retains reads with Ns in CIGAR strings, use <code>-U ALLOW_N_CIGAR_READS</code>.</li>
</ul>
<p>Subsetting reads corresponding to a genomic interval using PrintReads requires reads that are aligned to a reference genome, coordinate-sorted and indexed. Place the <code>.bai</code> index in the same directory as the <code>.bam</code> file.</p>
<pre><code class="pre_md">java -Xmx8G -jar /path/GenomeAnalysisTK.jar \
    -T PrintReads \ 
    -R /path/human_g1k_v37_decoy.fasta \ #reference fasta
    -L 10:91000000-92000000 \ #desired genomic interval chr:start-end
    -I 6517_2Mbp_input.bam \ #input
    -o 6517_1Mbp_output.bam </code class="pre_md"></pre>
<p>This creates a subset of reads from the input file, <code>6517_2Mbp_input.bam</code>, that align to the interval defined by the <code>-L</code> option, here a 1 Mbp region on chromosome 10. The tool creates two new files, <code>6517_1Mbp_output.bam</code> and corresponding index <code>6517_1Mbp_output.bai</code>. </p>
<ul>
<li>For paired reads, the tool does not account for reads whose mate aligns outside of the defined interval. To filter these lost mate reads, use RevertSam's <code>SANITIZE</code> option.</li>
</ul>
<p>To process large files, also designate a temporary directory. </p>
<pre><code class="pre_md">    TMP_DIR=/path/shlee #sets environmental variable for temporary directory</code class="pre_md"></pre>
<hr />