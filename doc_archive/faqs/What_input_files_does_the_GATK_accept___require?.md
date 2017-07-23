## What input files does the GATK accept / require?

http://gatkforums.broadinstitute.org/gatk/discussion/1204/what-input-files-does-the-gatk-accept-require

<p>All analyses done with the GATK typically involve several (though not necessarily all) of the following inputs:</p>
<ul>
<li>Reference genome sequence</li>
<li>Sequencing reads</li>
<li>Intervals of interest</li>
<li>Reference-ordered data</li>
</ul>
<p>This article describes the corresponding file formats that are acceptable for use with the GATK.</p>
<hr />
<h3>1. Reference Genome Sequence</h3>
<p>The GATK requires the reference sequence in a single reference sequence in FASTA format, with all contigs in the same file. The GATK requires strict adherence to the FASTA standard. All the standard IUPAC bases are accepted, but keep in mind that non-standard bases (i.e. other than ACGT, such as W for example) will be ignored (i.e. those positions in the genome will be skipped). </p>
<p><strong>Some users have reported having issues with reference files that have been stored or modified on Windows filesystems. The issues manifest as &quot;10&quot; characters (corresponding to encoded newlines) inserted in the sequence, which cause the GATK to quit with an error. If you encounter this issue, you will need to re-download a valid master copy of the reference file, or clean it up yourself.</strong> </p>
<p>Gzipped fasta files will not work with the GATK, so please make sure to unzip them first. Please see <a href="http://www.broadinstitute.org/gatk/guide/article?id=1601">this article</a> for more information on preparing FASTA reference sequences for use with the GATK.</p>
<h4>Important note about human genome reference versions</h4>
<p>If you are using human data, your reads must be aligned to one of the official b3x (e.g. b36, b37) or hg1x (e.g. hg18, hg19) references. The names and order of the contigs in the reference you used must exactly match that of one of the official references canonical orderings. These are defined by historical karotyping of largest to smallest chromosomes, followed by the X, Y, and MT for the b3x references; the order is thus 1, 2, 3, ..., 10, 11, 12, ... 20, 21, 22, X, Y, MT. The hg1x references differ in that the chromosome names are prefixed with &quot;chr&quot; and chrM appears first instead of last. The GATK will detect misordered contigs (for example, lexicographically sorted) and throw an error. This draconian approach, though unnecessary technically, ensures that all supplementary data provided with the GATK works correctly. You can use ReorderSam to fix a BAM file aligned to a missorted reference sequence.</p>
<p><strong>Our Best Practice recommendation is that you use a standard GATK reference from the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1213">GATK resource bundle</a>.</strong></p>
<hr />
<h3>2. Sequencing Reads</h3>
<p>The only input format for sequence reads that the GATK itself supports is the [Sequence Alignment/Map (SAM)] format. See [SAM/BAM] for more details on the SAM/BAM format as well as <a href="http://samtools.sourceforge.net/">Samtools</a> and <a href="http://picard.sourceforge.net/">Picard</a>, two complementary sets of utilities for working with SAM/BAM files.</p>
<p>If you don't find the information you need in this section, please see our <a href="http://www.broadinstitute.org/gatk/guide/article?id=1317">FAQs on BAM files</a>.</p>
<p>If you are starting out your pipeline with raw reads (typically in FASTQ format) you'll need to make sure that when you map those reads to the reference and produce a BAM file, the resulting BAM file is fully compliant with the GATK requirements. See the Best Practices documentation for detailed instructions on how to do this. </p>
<p>In addition to being in SAM format, we require the following additional constraints in order to use your file with the GATK:</p>
<ul>
<li>The file must be binary (with <code>.bam</code> file extension).</li>
<li>The file must be indexed.</li>
<li>The file must be sorted in coordinate order with respect to the reference (i.e. the contig ordering in your bam must exactly match that of the reference you are using).</li>
<li>The file must have a proper bam header with read groups. Each read group must contain the platform (PL) and sample (SM) tags. For the platform value, we currently support 454, LS454, Illumina, Solid, ABI_Solid, and CG (all case-insensitive).</li>
<li>Each read in the file must be associated with exactly one read group.</li>
</ul>
<p>Below is an example well-formed SAM field header and fields (with @SQ dictionary truncated to show only the first two chromosomes for brevity): </p>
<pre><code class="pre_md">@HD     VN:1.0  GO:none SO:coordinate
@SQ     SN:1    LN:249250621    AS:NCBI37       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:1b22b98cdeb4a9304cb5d48026a85128
@SQ     SN:2    LN:243199373    AS:NCBI37       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:a0d9851da00400dec1098a9255ac712e
@RG     ID:ERR000162    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR000252    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR001684    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@RG     ID:ERR001685    PL:ILLUMINA     LB:g1k-sc-NA12776-CEU-1 PI:200  DS:SRP000031    SM:NA12776      CN:SC
@PG     ID:GATK TableRecalibration      VN:v2.2.16      CL:Covariates=[ReadGroupCovariate, QualityScoreCovariate, DinucCovariate, CycleCovariate], use_original_quals=true, defau 
t_read_group=DefaultReadGroup, default_platform=Illumina, force_read_group=null, force_platform=null, solid_recal_mode=SET_Q_ZERO, window_size_nqs=5, homopolymer_nback=7, except on_if_no_tile=false, pQ=5, maxQ=40, smoothing=137       UR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta    M5:b4eb71ee878d3706246b7c1dbef69299
@PG     ID:bwa  VN:0.5.5
ERR001685.4315085       16      1       9997    25      35M     *       0       0       CCGATCTCCCTAACCCTAACCCTAACCCTAACCCT     ?8:C7ACAABBCBAAB?CCAABBEBA@ACEBBB@?     XT:A:U  XN:i:4    X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  RG:Z:ERR001685  NM:i:6  MD:Z:0N0N0N0N1A0A28     OQ:Z:&gt;&gt;:&gt;2&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;?&gt;&gt;&gt;&gt;??&gt;???&gt;
ERR001689.1165834       117     1       9997    0       *       =       9997    0       CCGATCTAGGGTTAGGGTTAGGGTTAGGGTTAGGG     &gt;7AA&lt;@@C?@?B?B??&gt;9?B??&gt;A?B???BAB??@     RG:Z:ERR001689    OQ:Z:&gt;:&lt;&lt;8&lt;&lt;&lt;&gt;&lt;&lt;&gt;&lt;&gt;&lt;&lt;&gt;7&lt;&gt;&gt;&gt;?&gt;&gt;??&gt;???????
ERR001689.1165834       185     1       9997    25      35M     =       9997    0       CCGATCTCCCTAACCCTAACCCTAACCCTAACCCT     758A:?&gt;&gt;8?=@@&gt;&gt;?;4&lt;&gt;=??@@==??@?==?8     XT:A:U  XN:i:4    SM:i:25 AM:i:0  X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  RG:Z:ERR001689  NM:i:6  MD:Z:0N0N0N0N1A0A28     OQ:Z:;74&gt;7&gt;&lt;&gt;&lt;&gt;&lt;&gt;&gt;&gt;&gt;&gt;&lt;:&lt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;&gt;
ERR001688.2681347       117     1       9998    0       *       =       9998    0       CGATCTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG     5@BA@A6B???A?B??&gt;B@B??&gt;B@B??&gt;BAB???     RG:Z:ERR001688    OQ:Z:=&gt;&gt;&gt;&gt;&lt;4&gt;&lt;&lt;?&gt;&lt;??????????????????????       </code class="pre_md"></pre>
<h4>Note about fixing BAM files with alternative sortings</h4>
<p>The GATK requires that the BAM file be sorted in the same order as the reference. Unfortunately, many BAM files have headers that are sorted in some other order -- lexicographical order is a common alternative. To resort the BAM file please use <a href="http://picard.sourceforge.net/command-line-overview.shtml#ReorderSam">ReorderSam</a>.   </p>
<hr />
<h3>3. Intervals of interest</h3>
<p>The GATK accept interval files for processing subsets of the genome in several different formats.  Please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1319">FAQs on interval lists</a> for details.</p>
<hr />
<h3>4. Reference Ordered Data (ROD) file formats</h3>
<p>The GATK can associate arbitrary reference ordered data (ROD) files with named tracks for all tools. Some tools require specific ROD data files for processing, and developers are free to write tools that access arbitrary data sets using the ROD interface. The general ROD system has the following syntax:</p>
<pre><code class="pre_md">-argumentName:name,type file</code class="pre_md"></pre>
<p>Where <code>name</code> is the name in the GATK tool (like &quot;eval&quot; in VariantEval), <code>type</code> is the type of the file, such as VCF or dbSNP, and <code>file</code> is the path to the file containing the ROD data.</p>
<p>The GATK supports several common file formats for reading ROD data:</p>
<ul>
<li><a href="http://www.1000genomes.org/wiki/analysis/variant-call-format/">VCF</a> : VCF type, the recommended format for representing variant loci and genotype calls. The GATK will only process valid VCF files; <a href="http://vcftools.sourceforge.net/">VCFTools</a> provides the official VCF validator. See <a href="http://vcftools.sourceforge.net/VCF-poster.pdf">here</a> for a useful poster detailing the VCF specification.</li>
<li>UCSC formated dbSNP : dbSNP type, UCSC dbSNP database output</li>
<li>BED : BED type, a general purpose format for representing genomic interval data, useful for masks and other interval outputs. <strong>Please note that the bed format is 0-based while most other formats are 1-based.</strong></li>
</ul>
<p><strong>Note that we no longer support the PED format. See <a href="http://atgu.mgh.harvard.edu/plinkseq/output.shtml">here</a> for converting .ped files to VCF.</strong></p>
<p>If you need additional information on VCF files, please see our FAQs on VCF files <a href="http://www.broadinstitute.org/gatk/guide/article?id=1318">here</a> and <a href="http://www.broadinstitute.org/gatk/guide/article?id=1268">here</a>.</p>