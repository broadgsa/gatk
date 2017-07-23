## (howto) Perform local realignment around indels

http://gatkforums.broadinstitute.org/gatk/discussion/7156/howto-perform-local-realignment-around-indels

<p><a name="top"></a>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/38/9f00ba2676151197dcf7a3ad96fe8f.png" align="right" height="150" style="margin:5px 0px 5px 10px"/>This tutorial replaces <a href="http://gatkforums.broadinstitute.org/discussion/2800/#top">Tutorial#2800</a> and applies to data types within the scope of the <a href="https://www.broadinstitute.org/gatk/guide/best-practices.php">GATK Best Practices</a> variant discovery workflow. </p>
<p>We provide example data and example commands for performing local realignment around small insertions and deletions (indels) against a reference. The resulting BAM reduces false positive SNPs and represents indels parsimoniously. First we use <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php">RealignerTargetCreator</a> to identify and create a target intervals list (<strong>step 1</strong>). Then we perform local realignment for the target intervals using <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php">IndelRealigner</a> (<strong>step 2</strong>). </p>
<hr />
<h4>Jump to a section</h4>
<ol>
<li><a href="#section0">Introduction</a></li>
<li><a href="#section1">Create target intervals list using RealignerTargetCreator</a></li>
<li><a href="#section2">Realign reads using IndelRealigner</a></li>
<li><a href="#section3">Some additional considerations</a> </li>
<li><a href="#section4">Related resources</a></li>
</ol>
<hr />
<p><a name="section0"></a></p>
<h2>1. Introduction and tutorial materials</h2>
<h4>Why do indel realignment?</h4>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/0f/721c2608b821393e108d8b0fd820ef.png" align="right" height="200" style="margin:5px 10px 5px 0px"/>Local realignment around indels allows us to correct mapping errors made by genome aligners and make read alignments more consistent in regions that contain indels.  </p>
<p>Genome aligners can only consider each read independently, and the scoring strategies they use to align reads relative to the reference limit their ability to align reads well in the presence of indels. Depending on the variant event and its relative location within a read, the aligner may favor alignments with mismatches or soft-clips instead of opening a gap in either the read or the reference sequence. In addition, the aligner's scoring scheme may use arbitrary tie-breaking, leading to different, non-parsimonious representations of the event in different reads.</p>
<p>In contrast, local realignment considers all reads spanning a given position. This makes it possible to achieve a high-scoring consensus that supports the presence of an indel event. It also produces a more parsimonious representation of the data in the region . </p>
<p>This two-step indel realignment process first identifies such regions where alignments may potentially be improved, then realigns the reads in these regions using a consensus model that takes all reads in the alignment context together.   </p>
<h4>Tools involved</h4>
<ul>
<li><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php">RealignerTargetCreator</a></li>
<li><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php">IndelRealigner</a></li>
</ul>
<h4>Prerequisites</h4>
<ul>
<li>Installed GATK tools</li>
<li>Coordinate-sorted and indexed BAM alignment data </li>
<li>Reference sequence, index and dictionary</li>
<li>An optional VCF file representing population variants, subset for indels</li>
</ul>
<h4>Download example data</h4>
<ul>
<li>To download the reference, open <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/">ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/</a> in your browser. Leave the password field blank. Download the following three files (~860 MB) to the same folder: <code>human_g1k_v37_decoy.fasta.gz</code>, <code>.fasta.fai.gz</code>, and <code>.dict.gz</code>. This same reference is available to load in IGV. </li>
<li>
<p>Click <a href="https://drive.google.com/open?id=0BzI1CyccGsZiQnBZdURMQkFobFk">tutorial_7156.tar.gz</a> to download the tutorial data. The data is human paired 2x150 whole genome sequence reads originally aligning at ~30x depth of coverage. The sample is a PCR-free preparation of the NA12878 individual run on the HiSeq X platform. I took the reads aligning to a one Mbp genomic interval (10:96,000,000-97,000,000) and sanitized and realigned the reads (BWA-MEM -M) to the entire genome according to the workflow presented in <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/">Tutorial#6483</a> and marked duplicates using MarkDuplicates according to <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6747/">Tutorial#6747</a>. We expect the alignment to reveal a good proportion of indels given its long reads (~150 bp per read), high complexity (PCR-free whole genome data) and deep coverage depth (30x). </p>
<p>Tutorial download also contains a known indels VCF from <a href="http://www.1000genomes.org/announcements/global-reference-human-genetic-variation-2015-09-30/">Phase 3 of the 1000 Genomes Project</a> subset for indel-only records in the interval 10:96,000,000-97,000,000. These represent consensus common and low-frequency indels in the studied populations from multiple approaches. The individual represented by our snippet, NA12878, is part of the 1000 Genomes Project data. Because of the differences in technology and methods used by the Project versus our sample library, our library has potential to reveal additional variants. </p>
</li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section1"></a></p>
<h2>2. Create target intervals list using RealignerTargetCreator</h2>
<p>For simplicity, we use a single known indels VCF, included in the tutorial data. For recommended resources, see <a href="https://www.broadinstitute.org/gatk/guide/article?id=1247">Article#1247</a>. </p>
<p>In the command, <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php">RealignerTargetCreator</a> takes a coordinate-sorted and indexed BAM and a VCF of known indels and creates a target intervals file.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R human_g1k_v37_decoy.fasta \
    -L 10:96000000-97000000 \
    -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    -I 7156_snippet.bam \
    -o 7156_realignertargetcreator.intervals</code class="pre_md"></pre>
<p>In the resulting file, <code>7156_realignertargetcreator.intervals</code>, intervals represent sites of extant and potential indels. If sites are proximal, the tool represents them as a larger interval spanning the sites. </p>
<h4>Comments on specific parameters</h4>
<ul>
<li>We specify the BAM alignment file with <code>-I</code>.</li>
<li>We specify the known indels VCF file with <code>-known</code>. The known indels VCF contains indel records only.</li>
<li>
<p>Three input choices are technically feasible in creating a target intervals list: you may provide RealignerTargetCreator (i) one or more VCFs of known indels each passed in via <code>-known</code>, (ii) one or more alignment BAMs each passed in via <code>-I</code> or (iii) both. We recommend the last mode, and we use it in the example command. We use these same input files again in the realignment step.</p>
<p>The tool adds indel sites present in the known indels file and indel sites in the alignment CIGAR strings to the targets. Additionally, the tool considers the presence of mismatches and soft-clips, and adds regions that pass a concentration threshold to the target intervals. </p>
<p>If you create an intervals list using only the VCF, RealignerTargetCreator will add sites of indel only records even if SNPs are present in the file. If you create an intervals list using both alignment and known indels, the known indels VCF should contain only indels. See <a href="#section4">Related resources</a>. </p>
</li>
<li>We include <code>-L 10:96000000-97000000</code> in the command to limit processing time. Otherwise, the tool traverses the entire reference genome and intervals outside these coordinates may be added given our example <code>7156_snippet.bam</code> contains a small number of alignments outside this region.</li>
<li>The tool samples to a target coverage of 1,000 for regions with greater coverage. </li>
</ul>
<p><a name="targetintervalsfile"></a></p>
<h4>The target intervals file</h4>
<p>The first ten rows of  <code>7156_realignertargetcreator.intervals</code> are as follows. The file is a text-based one-column list with one interval per row in <a href="https://www.biostars.org/p/84686/">1-based</a> coordinates. Header and column label are absent. For an interval derived from a known indel, the start position refers to the corresponding known variant. For example, for the first interval, we can <code>zgrep -w 96000399 INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf</code> for details on the 22bp deletion annotated at position 96000399. </p>
<pre><code>10:96000399-96000421
10:96002035-96002036
10:96002573-96002577
10:96003556-96003558
10:96004176-96004177
10:96005264-96005304
10:96006455-96006461
10:96006871-96006872
10:96007627-96007628
10:96008204</code></pre>
<p>To view intervals on IGV, convert the list to <a href="https://www.biostars.org/p/84686/">0-based</a> BED format using the following <a href="https://en.wikipedia.org/wiki/AWK">AWK</a> command. The command saves a new text-based file with <code>.bed</code> extension where chromosome, start and end are tab-separated, and the start position is one less than that in the intervals list.</p>
<pre><code class="pre_md">awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' 7156_realignertargetcreator.intervals &gt; 7156_realignertargetcreator.bed</code class="pre_md"></pre>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section2"></a></p>
<h2>3. Realign reads using IndelRealigner</h2>
<p>In the following command, <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php">IndelRealigner</a> takes a coordinate-sorted and indexed BAM and a target intervals file  generated by RealignerTargetCreator. IndelRealigner then performs local realignment on reads coincident with the target intervals using consenses from indels present in the original alignment.  </p>
<pre><code class="pre_md">java -Xmx8G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R human_g1k_v37_decoy.fasta \
    -targetIntervals 7156_realignertargetcreator.intervals \
    -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \ 
    -I 7156_snippet.bam \
    -o 7156_snippet_indelrealigner.bam</code class="pre_md"></pre>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/07/5c266d45fec220139c065ed8d540fb.png" align="right" height="330" style="margin:5px 0px 5px 10px"/> The resulting coordinate-sorted and indexed BAM contains the same records as the original BAM but with changes to realigned records and their mates. Our tutorial's two IGV screenshots show realigned reads in two different loci. For simplicity, the screenshots show the subset of reads that realigned. For screenshots of full alignments for the same loci, see <a href="https://us.v-cdn.net/5019796/uploads/FileUpload/f5/a3821923247a4ac14ce7fc554fbab2.png">here</a> and <a href="https://us.v-cdn.net/5019796/uploads/FileUpload/fe/dd603f62d298bae1c45cd2a9b36f75.png">here</a>.</p>
<h4>Comments on specific parameters</h4>
<ul>
<li>The <code>-targetIntervals</code> file from RealignerTargetCreator, with extension <code>.intervals</code> or <code>.list</code>, is required. See <a href="#targetintervalsfile">section 1</a> for a description.</li>
<li>Specify each BAM alignment file with <code>-I</code>. IndelRealigner operates on all reads simultaneously in files you provide it jointly. </li>
<li>Specify each optional known indels VCF file with <code>-known</code>. </li>
<li>For joint processing, e.g. for tumor-normal pairs, generate one output file for each input by specifying <code>-nWayOut</code> instead of <code>-o</code>. </li>
<li>
<p>By default, and in this command, IndelRealigner applies the <code>USE_READS</code> consensus model. This is the consensus model we recommend because it balances accuracy and performance. To specify a different model, use the <code>-model</code> argument. The <code>KNOWNS_ONLY</code> consensus model constructs alternative alignments from the reference sequence by incorporating any known indels at the site, the <code>USE_READS</code> model from indels in reads spanning the site and the <code>USE_SW</code> model additionally from Smith-Waterman alignment of reads that do not perfectly match the reference sequence. </p>
<p>The <code>KNOWNS_ONLY</code> model can be sufficient for preparing data for base quality score recalibration. It can maximize performance at the expense of some accuracy. This is the case only given the known indels file represents common variants for your data. If you specify <code>-model KNOWNS_ONLY</code> but forget to provide a VCF, the command runs but the tool does not realign any reads.   </p>
</li>
<li>If you encounter out of memory errors, try these options. First, increase max java heap size from <code>-Xmx8G</code>. To find a system's default maximum heap size, type <code>java -XX:+PrintFlagsFinal -version</code>, and look for <code>MaxHeapSize</code>. If this does not help, and you are jointly processing data, then try running indel realignment iteratively on smaller subsets of data before processing them jointly.</li>
<li>IndelRealigner performs local realignment without downsampling. If the number of reads in an interval exceeds the 20,000 default threshold set by the <code>-maxReads</code> parameter, then the tool skips the region.</li>
<li>The tool has two read filters, BadCigarFilter and MalformedReadFilter. The tool processes reads flagged as duplicate. </li>
</ul>
<h4>Changes to alignment records</h4>
<p>For our example data,194 alignment records realign for ~89 sites. These records now have the <code>OC</code> tag to mark the original CIGAR string. We can use the <code>OC</code> tag to pull out realigned reads and instructions for this are in <a href="#section4">section 4</a>. The following screenshot shows an example pair of records before and after indel realignment. We note seven changes with asterisks, blue for before and red for after, for both the realigned read and for its mate.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/d1/f2cfc19c2be272eb8b5eb03367e830.png" />
<p>Changes to the example realigned record:</p>
<ul>
<li>MAPQ increases from 60 to 70. The tool increases each realigned record's MAPQ by ten.</li>
<li>The <a href="http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F">CIGAR string</a>, now <code>72M20I55M4S</code>, reflects the realignment containing a 20bp insertion.</li>
<li>The OC tag retains the original CIGAR string (OC:Z:110M2I22M1D13M4S) and replaces the MD tag that stored the string for mismatching positions.</li>
<li>The NM tag counts the realigned record's mismatches, and changes from 8 to 24.</li>
</ul>
<p>Changes to the realigned read's mate record:</p>
<ul>
<li>The MC tag updates the mate CIGAR string (to MC:Z:72M20I55M4S).</li>
<li>The MQ tag updates to the new mapping quality of the mate (to MQ:i:70).</li>
<li>The UQ tag updates to reflect the new Phred likelihood of the segment, from UQ:i:100 to UQ:i:68.</li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section3"></a></p>
<h2>3. Some additional considerations</h2>
<p>RealignerTargetCreator documentation has a <code>-maxInterval</code> cutoff to drop intervals from the list if they are too large. This is because increases in number of reads per interval quadratically increase the compute required to realign a region, and larger intervals tend to include more reads. By the same reasoning, increasing read depth, e.g. with additional alignment files, increases required compute. </p>
<p>Our tutorial's <code>INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf</code> contains 1168 indel-only records. The following are metrics on intervals created using the three available options.</p>
<pre><code>               #intervals    avg length     basepair coverage     
VCF only       1161           3.33           3,864         
BAM only        487          15.22           7,412          
VCF+BAM        1151          23.07          26,558         </code></pre>
<p>You can project the genomic coverage of the intervals as a function of the interval density (number of intervals per basepair) derived from varying the known indel density (number of indel records in the VCF). This in turn allows you to anticipate compute for indel realignment. The density of indel sites increases the interval length following a power law (y=ax^b). The constant (a) and the power (b) are different for intervals created with VCF only and with VCF+BAM. For our example data, these average interval lengths are well within the length of a read and minimally vary the reads per interval and thus the memory needed for indel realignment. </p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section4"></a></p>
<h2>4. Related resources</h2>
<ul>
<li>See the <a href="https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1">Best Practice Workflow</a> and click on the flowchart's <code>Realign Indels</code> icon for best practice recommendations and links including to a 14-minute video overview.</li>
<li>See <a href="https://www.broadinstitute.org/gatk/guide/article?id=1247">Article#1247</a> for guidance on using VCF(s) of known variant sites. </li>
<li>To subset realigned reads only into a valid BAM, as shown in the screenshots, use <code>samtools view 7088_snippet_indelrealigner.bam | grep 'OC' | cut -f1 | sort &gt; 7088_OC.txt</code> to create a list of readnames. Then, follow direction in blogpost <a href="https://www.broadinstitute.org/gatk/blog?id=7019">SAM flags down a boat</a> on how to create a valid BAM using FilterSamReads. </li>
<li>See <a href="https://www.broadinstitute.org/gatk/guide/tagged?tag=multithreading">discussion on multithreading</a> for options on speeding up these processes. The document titled <a href="http://gatkforums.broadinstitute.org/gatk/discussion/1975/">How can I use parallelism to make GATK tools run faster?</a> gives two charts: (i) the first table relates the three parallelism options to the major GATK tools and (ii) the second table provides recommended configurations for the tools. Briefly, RealignerTargetCreator runs faster with increasing <code>-nt</code> threads, while IndelRealigner shows diminishing returns for increases in scatter-gather threads provided by Queue. See blog <a href="https://www.broadinstitute.org/gatk/blog?id=7249">How long does it take to run the GATK Best Practices?</a> for a breakdown of the impact of threading and CPU utilization for Best Practice Workflow tools. </li>
<li>See <a href="http://www.nature.com/ng/journal/v43/n5/full/ng.806.html">DePristo et al's 2011 <em>Nature Genetics</em> technical report</a> for benchmarked effects of indel realignment as well as for the mathematics behind the algorithms.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6517">Tutorial#6517</a> for instructions on creating a snippet of reads corresponding to a genomic interval. For your research aims, you may find testing a small interval of your alignment and your choice VCF, while adjusting parameters, before committing to processing your full dataset, is time well-invested. </li>
<li>The tutorial's PCR-free 2x150 bp reads give enough depth of coverage (34.67 mean and 99.6% above 15) and library complexity to allow us the confidence to use aligner-generated indels in realignment. Check alignment coverage with <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php">DepthofCoverage</a> for WGS or <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_diagnostics_diagnosetargets_DiagnoseTargets.php">DiagnoseTargets</a> for WES.</li>
<li>See <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php">SelectVariants</a> to subset out indel calls using the <code>-selectType INDEL</code> option. Note this excludes indels that are part of mixed variant sites (see <a href="http://gatkforums.broadinstitute.org/gatk/discussion/3682/">FAQ</a>). Current solutions to including indels from mixed sites involves the use of JEXL expressions, as discussed <a href="http://gatkforums.broadinstitute.org/gatk/discussion/1255/what-are-jexl-expressions-and-how-can-i-use-them-with-the-gatk">here</a>. Current solutions to selecting variants based on population allelic frequency (AF), as we may desire to limit our known indels to those that are more common than rare for more efficient processing, are discussed in two forum posts (<a href="http://gatkforums.broadinstitute.org/gatk/discussion/6526/selectvariants-af-with-multiallelic-variants">1</a>,<a href="http://gatkforums.broadinstitute.org/gatk/discussion/comment/21120#Comment_21120">2</a>). </li>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/6491/">Tutorial#6491</a> for basic instructions on using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>. </li>
</ul>
<p><a name="bottom"></a></p>