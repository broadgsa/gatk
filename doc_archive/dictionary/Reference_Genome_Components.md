## Reference Genome Components

http://gatkforums.broadinstitute.org/gatk/discussion/7857/reference-genome-components

<h4>Document is in <code>BETA</code>. It may be incomplete and/or inaccurate. Post suggestions to the <em>Comments</em> section.</h4>
<hr />
<p><a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a3/19c0c82fd10f748e04201847c89d70.png" align="right" height="210" style="margin:0px 0px 5px 10px"/></a> This document defines several components of a reference genome. We use the human GRCh38/hg38 assembly to illustrate.</p>
<p>GRCh38/hg38 is the assembly of the human genome released December of 2013, that uses alternate or <strong>ALT</strong> contigs to represent common complex variation, including <a href="https://en.wikipedia.org/wiki/Human_leukocyte_antigen">HLA</a> loci. Alternate contigs are also present in past assemblies but not to the extent we see with GRCh38. Much of the improvements in GRCh38 are the result of other genome sequencing and analysis projects, including the <a href="http://www.1000genomes.org/">1000 Genomes Project</a>. </p>
<p>The ideogram is from the <em>Genome Reference Consortium</em> website and showcases GRCh38.p7. The zoomed region illustrates how regions in blue are full of Ns.  </p>
<p><strong>Analysis set</strong> reference genomes have special features to accommodate sequence read alignment. This type of genome reference can differ from the reference you use to browse the genome.</p>
<ul>
<li>For example, the <a href="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/">GRCh38 analysis set</a> <strong>hard-masks</strong>, i.e. replaces with Ns, a proportion of homologous centromeric and genomic <a href="https://en.wikipedia.org/wiki/Satellite_DNA">repeat arrays</a> (on chromosomes 5, 14, 19, 21, &amp; 22) and two PAR (pseudoautosomal) regions on chromosome Y. Confirm the set you are using by viewing a PAR region of the Y chromosome on IGV as shown in the figure below. The chrY location of PAR1 and PAR2 on GRCh38 are chrY:10,000-2,781,479 and chrY:56,887,902-57,217,415.
<a href="https://us.v-cdn.net/5019796/uploads/FileUpload/83/c5938ded241dd754b8e8c148467338.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/83/c5938ded241dd754b8e8c148467338.png" align="" height="150" style=""/></a>
The sequence in the reference set is a mix of uppercase and lowercase letters. The lowercase letters represent <strong>soft-masked</strong> sequence corresponding to repeats from <a href="http://www.repeatmasker.org/">RepeatMasker</a> and <a href="https://tandem.bu.edu/trf/trf.html">Tandem Repeats Finder</a>.  </li>
<li>The GRCh38 analysis sets also include a contig to siphon off reads corresponding to the Epstein-Barr virus sequence as well as <strong>decoy</strong> contigs. The EBV contig can help correct for artifacts stemming from immortalization of human blood lymphocytes with <a href="https://en.wikipedia.org/wiki/Epstein%E2%80%93Barr_virus#Transformation_of_B-lymphocytes">EBV transformation</a>, as well as capture endogenous EBV sequence as <a href="http://gbe.oxfordjournals.org/content/6/4/846.full">EBV naturally infects B cells</a> in ~90% of the world population. Heng Li provides the decoy contigs.</li>
</ul>
<hr />
<h2>Nomenclature: words to describe components of reference genomes</h2>
<ul>
<li>
<p>A <strong>contig</strong> is a contiguous sequence without gaps.</p>
</li>
<li>
<p><strong>Alternate contigs</strong>, <strong>alternate scaffolds</strong> or <strong>alternate loci</strong> allow for representation of diverging haplotypes. These regions are too complex for a single representation. Identify ALT contigs by their <code>_alt</code> suffix.</p>
<p>The GRCh38 ALT contigs total 109Mb in length and span 60Mb of the primary assembly. Alternate contig sequences can be novel to highly diverged or nearly identical to corresponding primary assembly sequence. Sequences that are highly diverged from the primary assembly only contribute a few million bases. Most subsequences of ALT contigs are fairly similar to the primary assembly. This means that if we align sequence reads to GRCh38+ALT blindly, then we obtain many multi-mapping reads with zero mapping quality. Since many GATK tools have a ZeroMappingQuality filter, we will then miss variants corresponding to such loci.</p>
</li>
<li>
<p><strong>Primary assembly</strong> refers to the collection of (i) assembled chromosomes, (ii) unlocalized and (iii) unplaced sequences. It represents a non-redundant haploid genome.</p>
<p>(i) <strong>Assembled chromosomes</strong> for hg38 are chromosomes 1–22 (<code>chr1</code>–<code>chr22</code>), X (<code>chrX</code>), Y (<code>chrY</code>) and Mitochondrial (<code>chrM</code>).
(ii) <strong>Unlocalized</strong> sequence are on a specific chromosome but with unknown order or orientation. Identify by <code>_random</code> suffix.
(iii) <strong>Unplaced</strong> sequence are on an unknown chromosome. Identify by <code>chrU_</code> prefix.</p>
</li>
<li>
<p><strong>PAR</strong> stands for <a href="https://en.wikipedia.org/wiki/Pseudoautosomal_region">pseudoautosomal region</a>. PAR regions in mammalian X and Y chromosomes allow for recombination between the sex chromosomes. Because the PAR sequences together create a diploid or <em>pseudo-autosomal</em> sequence region, the X and Y chromosome sequences are intentionally identical in the genome assembly. <em>Analysis set</em> genomes further hard-mask two of the Y chromosome PAR regions so as to allow mapping of reads solely to the X chromosome PAR regions. </p>
</li>
<li>
<p>Different <strong>assemblies</strong> shift coordinates for loci and are released infrequently. Hg19 and hg38 represent two different major assemblies. Comparing data from different assemblies requires lift-over tools that adjust genomic coordinates to match loci, at times imperfectly. In the special case of hg19 and GRCh37, the primary assembly coordinates are identical for loci but patch updates differ. Also, the naming conventions of the references differ, e.g. the use of chr1 versus 1 to indicate chromosome 1, such that these also require lift-over to compare data. GRCh38/hg38 unifies the assemblies and the naming conventions.</p>
</li>
<li>
<p><strong>Patches</strong> are regional fixes that are released periodically for a given assembly. GRCh38.p7 indicates the seventh patched minor release of GRCh38. <a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/patches.shtml">This NCBI page</a> explains in more detail. Patches add information to the assembly without disrupting the chromosome coordinates. Again, they improve representation without affecting chromosome coordinate stability. The two types of patches, fixed and novel, represent different types of sequence.</p>
<p>(i) <strong>Fix patches</strong> represent sequences that will replace primary assembly sequence in the next major assembly release. When interpreting data, fix patches should take precedence over the chromosomes.
(ii) <strong>Novel patches</strong> represent alternate loci. When interpreting data, treat novel patches as population sequence variants.</p>
</li>
</ul>
<hr />
<h2>The GATK perspective on reference genomes</h2>
<p>Within GATK documentation, <a href="https://software.broadinstitute.org/gatk/documentation/article?id=8017">Tutorial#8017</a> outlines how to map reads in an alternate contig aware manner and discusses some of the implications of mapping reads to reference genomes with alternate contigs.   </p>
<p>GATK tools allow for use of a genomic <a href="http://gatkforums.broadinstitute.org/gatk/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals">intervals list</a> that tells tools which regions of the genome the tools should act on. Judicious use of an intervals list, e.g. one that excludes regions of Ns and low complexity repeat regions in the genome, makes processes more efficient. This brings us to the next point.</p>
<h4>Specifying contigs with colons in their names, as occurs for new contigs in GRCh38, requires special handling for GATK versions prior to v3.6. Please use the following workaround.</h4>
<ul>
<li>For example, <code>HLA-A*01:01:01:01</code> is a new contig in GRCh38. The colons are a new feature of contig naming for GRCh38 from prior assemblies. This has implications for using the <code>-L</code> option of GATK as the option also uses the colon as a delimiter to distinguish between contig and genomic coordinates.</li>
<li>When defining coordinates of interest for a contig, e.g. positions 1-100 for chr1, we would use <code>-L chr1:1-100</code>. This also works for our HLA contig, e.g. <code>-L HLA-A*01:01:01:01:1-100</code>.</li>
<li>
<p>However, when passing in an entire contig, for contigs with colons in the name, you must add <code>:1+</code> to the end of the chromosome name as shown below. This ensures that portions of the contig name are appropriately identified as part of the contig name and not genomic coordinates.</p>
<pre><code class="pre_md"> -L HLA-A*01:01:01:01:1+</code class="pre_md"></pre>
</li>
</ul>
<h3>Viewing CRAM alignments on genome browsers</h3>
<p>Because CRAM compression depends on the alignment reference genome, tools that use CRAM files ensure correct decompression by comparing reference contig <a href="https://en.wikipedia.org/wiki/MD5">MD5 hashtag</a> values. These are sensitive to any changes in the sequence, e.g. masking with Ns. This can have implications for viewing alignments in genome browsers when there is a disjoint between the reference that is loaded in the browser and the reference that was used in alignment. If you are using a version of tools for which this is an issue, be sure to load the original analysis set reference genome to view the CRAM alignments.</p>
<h3>Should I switch to a newer reference?</h3>
<p>Yes you should. In addition to adding many alternate contigs, GRCh38 corrects thousands of SNPs and indels in the GRCh37 assembly that are absent in the population and are likely sequencing artifacts. It also includes synthetic centromeric sequence and updates non-nuclear genomic sequence.</p>
<p>The ability to recognize alternate haplotypes for loci is a drastic improvement that GRCh38 makes possible. Going forward, expanding genomics data will help identify variants for alternate haplotypes, improve existing and add additional alternate haplotypes and give us a better accounting of alternate haplotypes within populations. We are already seeing improvements and additions in the patch releases to reference genomes, e.g. the seven minor releases of GRCh38 available at the time of this writing.  </p>
<p>Note that variants produced by alternate haplotypes when they are represented on the primary assembly may or may not be present in data resources, e.g. dbSNP. This could have varying degrees of impact, including negligible, for any process that relies on known variant sites. Consider the impact this discrepant coverage in data resources may have for your research aims and weigh this against the impact of missing variants because their sequence context is unaccounted for in previous assemblies.</p>
<hr />
<h2>External resources</h2>
<ol>
<li><code>New 11/16/2016</code> For a brief history and discussion on challenges in using GRCh38, see the 2015 <em>Genome Biology</em> article <em>Extending reference assembly models</em> by Church et al. (DOI: <a href="https://dx.doi.org/10.1186/s13059-015-0587-3">10.1186/s13059-015-0587-3</a>).</li>
<li>For press releases highlighting improvements in GRCh38 from December 2013, see <a href="http://www.ncbi.nlm.nih.gov/news/12-23-2013-grch38-released/">http://www.ncbi.nlm.nih.gov/news/12-23-2013-grch38-released/</a> and <a href="http://genomeref.blogspot.co.uk/2013/12/announcing-grch38.html">http://genomeref.blogspot.co.uk/2013/12/announcing-grch38.html</a>. The latter post summarizes major improvements, including the correction of thousands of SNPs and indels in GRCh37 not seen in the population and the inclusion of synthetic centromeric sequence.</li>
<li>Recent releases of BWA, e.g. v0.7.15+, handle alt contig mapping and HLA typing. See the <a href="https://github.com/lh3/bwa">BWA repository</a> for information. See these pages for <a href="https://sourceforge.net/projects/bio-bwa/files/">download</a> and <a href="http://gatkforums.broadinstitute.org/wdl/discussion/2899">installation instructions</a>.</li>
<li>The <a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/">Genome Reference Consortium (GRC)</a> provides human, mouse, zebrafish and chicken sequences, and <a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/">this particular webpage</a> gives an overview of GRCh38. Namely, an interactive chromosome ideogram marks regions with corresponding alternate loci, regions with fix patches and regions containing novel patches. For additional assembly terminology, see <a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/definitions.shtml"><a href="http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/definitions.shtml">http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/definitions.shtml</a></a>.</li>
<li>
<p>The <a href="http://genome.ucsc.edu/cgi-bin/hgGateway?clade=mammal&amp;org=Human&amp;db=hg38">UCSC Genome Browser</a> allows browsing and download of genomes, including <em>analysis sets</em>, from many different species. For more details on the difference between GRCh38 reference and analysis sets, see <code>ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/README.txt</code> and <code>ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/README.txt</code>, respectively. In addition, the site provides annotation files, e.g. <a href="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/">here</a> is the annotation database for GRCh38. Within this particular page, the file named <a href="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz">gap.txt.gz</a> catalogues the gapped regions of the assembly full of Ns. For our illustration above, the corresponding region in this file shows:</p>
<pre><code class="pre_md">    585    chr14    0    10000    1    N    10000    telomere    no
    1    chr14    10000    16000000    2    N    15990000    short_arm    no
    707    chr14    16022537    16022637    4    N    100    contig    no</code class="pre_md"></pre>
</li>
<li>The <a href="http://www.broadinstitute.org/igv/home">Integrative Genomics Viewer</a> is a desktop application for viewing genomics data including alignments. The tool accesses reference genomes you provide via file or URL or that it hosts over a server. The numerous hosted reference genomes include GRCh38. See <a href="http://www.broadinstitute.org/igv/Genomes">this page</a> for information on hosted reference genomes. For the most up-to-date list of hosted genomes, open IGV and go to <em>Genomes</em>&gt;<em>Load Genome From Server</em>. A menu lists genomes you can make available in the main genome dropdown menu. </li>
</ol>
<hr />