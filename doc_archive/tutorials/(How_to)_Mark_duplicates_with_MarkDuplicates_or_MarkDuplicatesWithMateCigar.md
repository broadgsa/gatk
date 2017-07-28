## (How to) Mark duplicates with MarkDuplicates or MarkDuplicatesWithMateCigar

http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar

<p><a name="top"></a>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/fc/ae1f9a500dd04a157a46b653c8f1e2.png" height="180"align="right" border="9"/> This tutorial updates <a href="http://gatkforums.broadinstitute.org/discussion/2799/#top">Tutorial#2799</a>.</p>
<p>Here we discuss two tools, <a href="http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates">MarkDuplicates</a> and <a href="http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicatesWithMateCigar">MarkDuplicatesWithMateCigar</a>, that flag duplicates. We provide example data and example commands for you to follow along the tutorial (<strong>section 1</strong>) and include tips in estimating library complexity for PCR-free samples and patterned flow cell technologies. In <strong>section 2</strong>, we point out special memory considerations for these tools. In <strong>section 3</strong>, we highlight the similarities and differences between the two tools. Finally, we get into some details that may be of interest to some that includes comments on the metrics file (<strong>section 4</strong>).</p>
<blockquote>
<p>To mark duplicates in RNA-Seq data, use MarkDuplicates. Reasons are explained in <a href="#section2">section 2</a> and <a href="#section3">section 3</a>. And if you are considering using MarkDuplicatesWithMateCigar for your DNA data, be sure insert lengths are short and you have a low percentage of split or multi-mapping records.</p>
</blockquote>
<p>Obviously, expect more duplicates for samples prepared with PCR than for PCR-free preparations. Duplicates arise from various sources, including within the sequencing run. As such, even PCR-free data can give rise to duplicates, albeit at low rates, as illustrated here with our example data.</p>
<h3>Which tool should I use, MarkDuplicates or MarkDuplicatesWithMateCigar? <code>new section 5/25/2016</code></h3>
<p>The Best Practices so far recommends MarkDuplicates. However, as always, consider your research goals.</p>
<p>If your research uses paired end reads and pre-processing that generates missing mates, for example by application of an intervals list or by removal of reference contigs after the initial alignment, and you wish to flag duplicates for these remaining singletons, then MarkDuplicatesWithMateCigar will flag these for you at the insert level using the mate cigar (MC) tag. MarkDuplicates skips these singletons from consideration. </p>
<p>If the qualities by which the representative insert in a duplicate set is selected is important to your analyses, then note that MarkDuplicatesWithMateCigar is limited to prioritizing by the total mapped length of a pair, while MarkDuplicates can use this OR the default sum of base qualities of a pair. </p>
<p>If you are still unsure which tool is appropriate, then consider maximizing comparability to previous analyses. The Broad Genomics Platform has used only MarkDuplicates in their production pipelines. MarkDuplicatesWithMateCigar is a newer tool that has yet to gain traction. </p>
<p>This tutorial compares the two tools to dispel the circulating notion that the outcomes from the two tools are equivalent and to provide details helpful to researchers in optimizing their analyses.</p>
<p>We welcome feedback. Share your suggestions in the <a href="#bottom">Comment section</a> at the bottom of this page.</p>
<hr />
<h4>Jump to a section</h4>
<ol>
<li><a href="#section1">Commands for MarkDuplicates and MarkDuplicatesWithMateCigar</a></li>
<li><a href="#section2">Slow or <em>out of memory</em> error? Special memory considerations for duplicate marking tools</a> </li>
<li><a href="#section3">Conceptual overview of duplicate flagging</a> </li>
<li><a href="#section4">Details of interest to some</a></li>
</ol>
<hr />
<h4>Tools involved</h4>
<ul>
<li><a href="http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates">MarkDuplicates</a></li>
<li><a href="http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicatesWithMateCigar">MarkDuplicatesWithMateCigar</a></li>
</ul>
<h4>Prerequisites</h4>
<ul>
<li>Installed Picard tools</li>
<li>Coordinate-sorted and indexed BAM alignment data. Secondary/supplementary alignments are flagged appropriately (256 and 2048 flags) and additionally with the mate unmapped (8) flag. See the <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/#step3C">MergeBamAlignment section (3C)</a> of <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483">Tutorial#6483</a> for a description of how MergeBamAlignment ensures such flagging. <code>**Revision as of 5/17/2016:**</code> I wrote this tutorial at a time when the input could only be an indexed and coordinate-sorted BAM. Recently, the tools added a feature to accept queryname-sorted inputs that in turn activates additional features. The additional features that providing a queryname-sorted BAM activates will give DIFFERENT duplicate flagging results. So for the tutorial's observations to apply, use coordinate-sorted data. </li>
<li>For MarkDuplicatesWithMateCigar, pre-computed Mate CIGAR (MC) tags. Data produced according to <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/">Tutorial#6483</a> will have the MC tags added by MergeBamAlignment. Alternatively, see tools <a href="http://broadinstitute.github.io/picard/command-line-overview.html#RevertOriginalBaseQualitiesAndAddMateCigar">RevertOriginalBaseQualitiesAndAddMateCigar</a> and <a href="http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation">FixMateInformation</a>.</li>
<li>Appropriately assigned Read Group (RG) information. Read Group library (RGLB) information is factored for molecular duplicate detection. Optical duplicates are limited to those from the same RGID. </li>
</ul>
<h4>Download example data</h4>
<ul>
<li>Use the <a href="http://gatkforums.broadinstitute.org/discussion/4610/">advanced tutorial bundle</a>'s human_g1k_v37_decoy.fasta as reference. This same reference is available to load in IGV. </li>
<li><a href="https://drive.google.com/open?id=0BzI1CyccGsZiWURLdUdfRjVQazg">tutorial_6747.tar.gz</a> data contain human paired 2x150 whole genome sequence reads originally aligning at ~30x depth of coverage. The sample is a PCR-free preparation of the NA12878 individual run on the HiSeq X platform. This machine type, along with HiSeq 4000, has the newer patterned flow cell that differs from the typical non-patterned flow cell. I took the reads aligning to a one Mbp genomic interval (10:96,000,000-97,000,000) and sanitized and realigned the reads (BWA-MEM -M) to the entire genome according to the workflow presented in <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/">Tutorial#6483</a> to produce <code>snippet.bam</code>. The data has (i) no supplementary records; (ii) secondary records flagged with the 256 flag <em>and</em> the mate-unmapped (8) flag; and (iii) unmapped records (4 flag) with mapped mates (mates have 8 flag), zero MAPQ (column 5) and asterisks for CIGAR (column 6). The notation allows read pairs where one mate maps and the other does not to sort and remain together when we apply genomic intervals such as in the generation of the snippet.</li>
</ul>
<h4>Related resources</h4>
<ul>
<li>See <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics">DuplicationMetrics</a> for descriptions of each metric.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/">Tutorial#6483</a> for instructions on how to efficiently map and clean up short read sequence data. You can use the resulting files directly in this tutorial.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/firecloud/discussion/3059/lane-library-sample-and-cohort-what-do-they-mean-and-why-are-they-important">an overview of lane, library, sample and cohort</a> and <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6199/picard-mark-duplicates-handling-of-library-information">this forum discussion of how MarkDuplicates handles library information</a>.</li>
<li>See <a href="https://broadinstitute.github.io/picard/explain-flags.html">SAM flags</a> to interpret SAM flag values.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6329">dictionary entry on Illumina Chastity filter</a> for a link to a document comparing patterned and non-patterned flow cells. </li>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/2909/">this tutorial</a> to coordinate-sort and index a BAM.</li>
<li>See <a href="http://gatkforums.broadinstitute.org/discussion/6491/">this tutorial</a> for basic instructions on using the <a href="http://www.broadinstitute.org/igv/">Integrative Genomics Viewer (IGV)</a>. </li>
</ul>
<hr />
<p><a name="section1"></a></p>
<h2>1. Commands for MarkDuplicates and MarkDuplicatesWithMateCigar</h2>
<p>The following commands take a coordinate-sorted and indexed BAM and return (i) a BAM with the same records in coordinate order and with duplicates marked by the 1024 flag, (ii) a duplication metrics file, and (iii) an optional matching BAI index. </p>
<p>For a given file with all MC (mate CIGAR) tags accounted for:</p>
<ul>
<li>and where all mates are accounted for, each tool--MarkDuplicates and MarkDuplicatesWithMateCigar--examines the same duplicate sets but prioritize which inserts get marked duplicate differently. This situation is represented by our <code>snippet</code> example data.</li>
<li>but containing missing mates records, MarkDuplicates ignores the records, while MarkDuplicatesWithMateCigar still considers them for duplicate marking using the MC tag for mate information. Again, the duplicate scoring methods differ for each tool. </li>
</ul>
<p>Use the following commands to flag duplicates for <code>6747_snippet.bam</code>. These commands produce qualitatively different data.</p>
<p><strong>Score duplicate sets based on the sum of base qualities using MarkDuplicates:</strong></p>
<pre><code class="pre_md">java -Xmx32G -jar picard.jar MarkDuplicates \
INPUT=6747_snippet.bam \ #specify multiple times to merge 
OUTPUT=6747_snippet_markduplicates.bam \
METRICS_FILE=6747_snippet_markduplicates_metrics.txt \ 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \ #changed from default of 100
CREATE_INDEX=true \ #optional
TMP_DIR=/tmp</code class="pre_md"></pre>
<p><strong>Score duplicate sets based on total mapped reference length using MarkDuplicatesWithMateCigar:</strong></p>
<pre><code class="pre_md">java -Xmx32G -jar picard.jar MarkDuplicatesWithMateCigar \
INPUT=6747_snippet.bam \ #specify multiple times to merge
OUTPUT=6747_snippet_markduplicateswithmatecigar.bam \
METRICS_FILE=6747_snippet_markduplicateswithmatecigar_metrics.txt \ 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \ #changed from default of 100
CREATE_INDEX=true \ #optional
TMP_DIR=/tmp</code class="pre_md"></pre>
<h3>Comments on select parameters</h3>
<ul>
<li><code>**Revision as of 5/17/2016:**</code> The example input <code>6747_snippet.bam</code> is coordinate-sorted and indexed. Recently, the tools added a feature to accept queryname-sorted inputs that in turn by default activates additional features that will give DIFFERENT duplicate flagging results than outlined in this tutorial. Namely, if you provide MarkDuplicates a queryname-sorted BAM, then if a primary alignment is marked as duplicate, then the tool will also flag its (i) unmapped mate, (ii) secondary and/or (iii) supplementary alignment record(s) as duplicate.</li>
<li>Each tool has a distinct default <code>DUPLICATE_SCORING_STRATEGY</code>. For MarkDuplicatesWithMateCigar it is TOTAL_MAPPED_REFERENCE_LENGTH and this is the <em>only</em> scoring strategy available. For MarkDuplicates you can switch the <code>DUPLICATE_SCORING_STRATEGY</code> between the default SUM_OF_BASE_QUALITIES and TOTAL_MAPPED_REFERENCE_LENGTH. Both scoring strategies use <em>information</em> pertaining to both mates in a pair, but in the case of MarkDuplicatesWithMateCigar the information for the mate comes from the read's MC tag and not from the actual mate. </li>
<li>To <strong>merge multiple files into a single output</strong>, e.g. when aggregating a sample from across lanes, specify the <code>INPUT</code> parameter for each file. The tools merge the read records from the multiple files into the single output file. The tools marks duplicates for the entire library (RGLB) and accounts for optical duplicates per RGID. <code>INPUT</code> files must be coordinate sorted and indexed.</li>
<li>The Broad's production workflow increases  <code>OPTICAL_DUPLICATE_PIXEL_DISTANCE</code> to 2500, to better estimate library complexity. The default setting for this parameter is 100. Changing this parameter does not alter duplicate marking. It only changes the count for optical duplicates and the library complexity estimate in the metrics file in that whatever is counted as an optical duplicate does not factor towards library complexity. The increase has to do with the fact that our example data was sequenced in a patterned flow cell of a HiSeq X machine. Both HiSeq X and HiSeq 4000 technologies decrease pixel unit area by 10-fold and so the equivalent pixel distance in non-patterned flow cells is 250. You may ask why are we still counting optical duplicates for patterned flow cells that by design should have no optical duplicates. We are hijacking this feature of the tools to account for other types of duplicates arising from the sequencer. Sequencer duplicates are not limited to optical duplicates and should be differentiated from PCR duplicates for more accurate library complexity estimates. </li>
<li>By default the tools flag duplicates and retain them in the output file. <strong>To remove the duplicate records</strong> from the resulting file, set the <code>REMOVE_DUPLICATES</code> parameter to true. However, given you can set GATK tools to include duplicates in analyses by adding <code>-drf DuplicateRead</code> to commands, a better option for value-added storage efficiency is to retain the resulting marked file over the input file.</li>
<li>To <strong>optionally create a <code>.bai</code> index</strong>, add and set the <code>CREATE_INDEX</code> parameter to true. </li>
</ul>
<p>For <code>snippet</code>, the duplication metrics are identical whether marked by MarkDuplicates or MarkDuplicatesWithMateCigar. We have 13.4008% duplication, with 255 <em>unpaired read duplicates</em> and 18,254 <em>paired read duplicates</em>. However, as the screenshot at the top of this page illustrates, and as <a href="#section4">section 4</a> explains, the data qualitatively differ. </p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section2"></a></p>
<h2>2. Slow or <em>out of memory</em> error? Special memory considerations for duplicate marking tools</h2>
<p>The seemingly simple task of marking duplicates is one of the most memory hungry processes, especially for paired end reads. Both tools are compute-intensive and require upping memory compared to other processes. </p>
<p>Because of the single-pass nature of MarkDuplicatesWithMateCigar, for a given file its memory requirements can be greater than for MarkDuplicates. What this means is that MarkDuplicatesWithMateCigar streams the duplicate marking routine in a manner that allows for <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6483/#step3D">piping</a>. Due to these memory constraints for MarkDuplicatesWithMateCigar, we recommend MarkDuplicates for alignments that have large reference skips, e.g. spliced RNA alignments. </p>
<p>For large files, (1) use the Java <code>-Xmx</code> setting and (2) set the environmental variable <code>TMP_DIR</code> for a temporary directory. These options allow the tool to run without slowing down as well as run without causing an <em>out of memory</em> error. For the purposes of this tutorial, commands are given as if the example data is a large file, which we know it is not. </p>
<pre><code class="pre_md">    java -Xmx32G -jar picard.jar MarkDuplicates \
    ... \
    TMP_DIR=/tmp </code class="pre_md"></pre>
<p>These options can be omitted for small files such as the example data and the equivalent command is as follows.</p>
<pre><code class="pre_md">    java -jar picard.jar MarkDuplicates ...   </code class="pre_md"></pre>
<h3>Set the java maxheapsize, specified by the <code>-Xmx#G</code> option, to the maximum your system allows.</h3>
<p>The high memory cost, especially for MarkDuplicatesWithMateCigar, is in part because the tool systematically traverses genomic coordinate intervals for inserts in question, and for every read it marks as a duplicate it must keep track of the mate, which may or may not map nearby, so that reads are marked as pairs with each record emitted in its coordinate turn. In the meanwhile, this information is held in memory, which is the first choice for faster processing, until the memory limit is reached, at which point memory spills to disk. We set this limit high to minimize instances of memory spilling to disk.</p>
<p>In the example command, the <code>-Xmx32G</code> Java option caps the maximum heap size, or memory usage, to 32 gigabytes, which is the limit on the server I use. This is in contrast to the 8G setting I use for other processes on the same sample data--a 75G BAM file. To find a system's default maximum heap size, type <code>java -XX:+PrintFlagsFinal -version</code>, and look for <code>MaxHeapSize</code>. </p>
<h3>Set an additional temporary directory with the <code>TMP_DIR</code> parameter for memory spillage.</h3>
<p>When the tool hits the memory limit, memory spills to disk. This causes data to traverse in and out of the processor's I/O device, slowing the process down. Disk is a location you specify with the <code>TMP_DIR</code> parameter. If you work on a server separate from where you read and write files to, setting TMP_DIR to the server's local temporary directory (typically <code>/tmp</code>) can reduce processing time compared to setting it to the storage disk. This is because the tool then additionally avoids traversing the network file system when spilling memory. Be sure the TMP_DIR location you specify provides enough storage space. Use <code>df -h</code> to see how much is available. </p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section3"></a></p>
<h2>3. Conceptual overview of duplicate flagging</h2>
<p><strong>The aim of duplicate marking</strong> is to flag all but one of a duplicate set as duplicates and to use duplicate metrics to estimate library complexity. Duplicates have a higher probability of being non-independent measurements from the exact same template DNA. Duplicate inserts are marked by the 0x400 bit (1024 flag) in the second column of a SAM record, for each mate of a pair. This allows downstream GATK tools to exclude duplicates from analyses (most do this by default). Certain duplicates, i.e. PCR and sequencer duplicates, violate assumptions of variant calling and also potentially amplify errors. Removing these, even at the cost of removing serendipitous biological duplicates, allows us to be conservative in calculating the confidence of variants.  </p>
<blockquote>
<p>GATK tools allow you to <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_CommandLineGATK.php">disable the duplicate read filter</a> with <code>-drf DuplicateRead</code> so you can include duplicates in analyses. </p>
</blockquote>
<p>For a whole genome DNA sample, <strong>duplicates arise from three sources</strong>: (i) in DNA shearing from distinct molecular templates identical in insert mapping, (ii) from PCR amplification of a template (PCR duplicates), and (iii) from sequencing, e.g. optical duplicates. The tools cannot distinguish between these types of duplicates with the exception of optical duplicates. In estimating library complexity, the latter two types of duplicates are undesirable and should each factor differently. </p>
<p><strong>When should we not care about duplicates?</strong> Given duplication metrics, we can make some judgement calls on the quality of our sample preparation and sequencer run. Of course, we may not expect a complex library if our samples are targeted amplicons. Also, we may expect minimal duplicates if our samples are PCR-free. Or it may be that because of the variation inherent in expression level data, e.g. RNA-Seq, duplicate marking becomes ritualistic. Unless you are certain of your edge case (amplicon sequencing, RNA-Seq allele-specific expression analysis, etc.) where duplicate marking adds minimal value, you should go ahead and mark duplicates. You may find yourself staring at an IGV session trying to visually calculate the strength of the evidence for a variant. We can pat ourselves on the back for having the forethought to systematically mark duplicates and turn on the IGV duplicate filter. </p>
<blockquote>
<p>The <strong>Broad's Genomics Platform uses MarkDuplicates twice for multiplexed samples</strong>. Duplicates are flagged first per sample per lane to estimate lane-level library complexity, and second to aggregate data per sample while marking all library duplicates. In the second pass, duplicate marking tools again assess all reads for duplicates and overwrite any prior flags.</p>
</blockquote>
<p>Our two duplicate flagging tools <strong>share common features but differ at the core</strong>. As the name implies, MarkDuplicatesWithMateCigar uses the MC (mate CIGAR) tag for mate alignment information. Unlike MarkDuplicates, it is a single-pass tool that requires pre-computed MC tags.</p>
<ul>
<li>For RNA-Seq data mapped against the genome, use MarkDuplicates. Specifically, MarkDuplicatesWithMateCigar will refuse to process data with large reference skips frequent in spliced RNA transcripts where the gaps are denoted with an <code>N</code> in the CIGAR string. </li>
<li>Both tools only consider primary mappings, even if mapped to different contigs, and ignore secondary/supplementary alignments (256 flag and 2048 flag) altogether.  Because of this, before flagging duplicates, be sure to mark primary alignments according to a strategy most suited to your experimental aims. See <a href="http://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a>'s <code>PRIMARY_ALIGNMENT_STRATEGY</code> parameter for strategies the tool considers for changing primary markings made by an aligner. </li>
<li>Both tools identify duplicate sets identically with the exception that MarkDuplicatesWithMateCigar additionally considers reads with missing mates. Missing mates occur for example when aligned reads are filtered using an interval list of genomic regions. This creates divorced reads whose mates aligned outside the targeted intervals.</li>
<li>Both tools identify duplicates as sets of read pairs that have the same unclipped alignment start and unclipped alignment end. The tools intelligently factor for discordant pair orientations given these start and end coordinates. Within a duplicate set, with the exception of optical duplicates, read pairs may have either pair orientation--F1R2 or F2R1. For optical duplicates, pairs in the set must have the same orientation. Why this is is explained in <a href="#section4">section 4</a>.</li>
<li>Both tools take into account clipped and gapped alignments and singly mapping reads (mate unmapped and not secondary/supplementary). </li>
<li>Each tool flags duplicates according to different priorities. MarkDuplicatesWithMateCigar prioritizes which pair to leave as the representative nondup based on the total mapped length of a pair while MarkDuplicates can prioritize based on the sum of base qualities of a pair (default) or the total mapped length of a pair. Duplicate <em>inserts</em> are marked at both ends. </li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="section4"></a></p>
<h2>4. Details of interest to some</h2>
<p>To reach a high target coverage depth, some fraction of sequenced reads will by stochastic means be duplicate reads. </p>
<p>Let us hope the truth of a variant never comes down to so few reads that duplicates should matter so. Keep in mind the better evidence for a variant is the presence of overlapping reads that contain the variant. Also, take estimated library complexity at face value--an estimate. </p>
<h3>Don't be duped by identical numbers. Data from the two tools <em>qualitatively differ</em>.</h3>
<p>First, let me reiterate that secondary and supplementary alignment records are skipped and never flagged as duplicate. </p>
<p>Given a file with no missing mates, each tool identifies the same duplicate sets from primary alignments only and therefore the <em>same number</em> of duplicates. To reiterate, the number of identical loci or duplicate sets and the records within each set are the same for each tool. However, each tool differs in how it decides which insert(s) within a set get flagged and thus which insert remains the representative <em>nondup</em>. Also, if there are ties, the tools may break them differently in that tie-breaking can depend on the sort order of the records in memory. </p>
<ul>
<li>MarkDuplicates by default prioritizes the sum of base qualities for both mates of a pair. The pair with the highest sum of base qualities remains as the nondup. </li>
<li>As a consequence of using the mate's CIGAR string (provided by the MC tag), MarkDuplicatesWithMateCigar can only prioritize the total mapped reference length, as provided by the CIGAR string, in scoring duplicates in a set. The pair with the longest mapping length remains as the nondup. </li>
<li>If there are ties after applying each scoring strategy, both tools break the ties down a chain of deterministic factors starting with read name.</li>
</ul>
<h3>Duplicate metrics in brief</h3>
<p>We can break down the metrics file into two parts: (1) a table of metrics that counts various categories of duplicates and gives the library complexity estimate, and (2) histogram values in two columns.</p>
<p>See <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics">DuplicationMetrics</a> <strong>for descriptions of each metric</strong>. For paired reads, duplicates are considered for the insert. For single end reads, duplicates are considered singly for the read, increasing the likelihood of being identified as a duplicate. Given the lack of insert-level information for these singly mapping reads, the insert metrics calculations exclude these. </p>
<p>The <strong>library complexity estimate</strong> only considers the duplicates that remain after subtracting out optical duplicates. For the math to derive estimated library size, see formula (1.2) in <a href="https://www.broadinstitute.org/gatk/media/docs/Samtools.pdf">Mathematical Notes on SAMtools Algorithms</a>. </p>
<p>The <strong>histogram values</strong> extrapolate the calculated library complexity to a saturation curve plotting the gains in complexity if you sequence additional aliquots of the same library. The first bin's value represents the current complexity. </p>
<h3>Pair orientation F1R2 is distinct from F2R1 for optical duplicates</h3>
<p>Here we refer you to a <a href="https://www.youtube.com/watch?v=womKfikWlxM">five minute video</a> illustrating what happens at the molecular level in a typical sequencing by synthesis run.</p>
<p>What I would like to highlight is that each strand of an insert has a chance to seed a different cluster. I will also point out, due to sequencing chemistry, F1 and R1 reads typically have better base qualities than F2 and R2 reads. </p>
<blockquote>
<p>Optical duplicate designation requires the same pair orientation.</p>
</blockquote>
<p>Let us work out the implications of this for a paired end, unstranded DNA library. During sequencing, within the flow cell, for a particular insert produced by sample preparation, the strands of the insert are separated and each strand has a chance to seed a different cluster. Let's say for InsertAB, ClusterA and ClusterB and for InsertCD, ClusterC and ClusterD. InsertAB and InsertCD are identical in sequence and length and map to the same loci. It is possible InsertAB and InsertCD are PCR duplicates and also possible they represent original inserts. Each strand is then sequenced in the forward and reverse to give four pieces of information in total for the given insert, e.g. ReadPairA and ReadPairB for InsertAB. The pair orientation of these two pairs are reversed--one cluster will give F1R2 and the other will give F2R1 pair orientation. Both read pairs map exactly to the same loci. Our duplicate marking tools consider ReadPairA and ReadPairB in the same duplicate set for regular duplicates but not for optical duplicates. Optical duplicates require identical pair orientation.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="bottom"></a></p>