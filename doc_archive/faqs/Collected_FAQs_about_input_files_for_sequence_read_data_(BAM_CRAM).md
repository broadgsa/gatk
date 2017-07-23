## Collected FAQs about input files for sequence read data (BAM/CRAM)

http://gatkforums.broadinstitute.org/gatk/discussion/1317/collected-faqs-about-input-files-for-sequence-read-data-bam-cram

<h3>1. What file formats do you support for sequence data input?</h3>
<p>The GATK supports the <a href="http://samtools.github.io/hts-specs/">BAM</a> format for reads, quality scores, alignments, and metadata (<em>e.g.</em> the lane of sequencing, center of origin, sample name, etc.). Starting with version 3.5, the <a href="http://samtools.github.io/hts-specs/">CRAM</a> format is supported as well. SAM format is not supported but can be easily converted with Picard tools. </p>
<hr />
<h3>2. How do I get my data into BAM format?</h3>
<p>The GATK doesn't have any tools for getting data into BAM format, but many other toolkits exist for this purpose. We recommend you look at <a href="http://broadinstitute.github.io/picard/">Picard</a> and <a href="http://samtools.sourceforge.net/">Samtools</a> for creating and manipulating BAM files. Also, many aligners are starting to emit BAM files directly. See <a href="http://bio-bwa.sourceforge.net/bwa.shtml">BWA</a> for one such aligner.</p>
<hr />
<h3>3. What are the formatting requirements for my BAM file(s)?</h3>
<p>All BAM/CRAM files must satisfy the following requirements:</p>
<ul>
<li>It must be aligned to one of the references described <a href="http://www.broadinstitute.org/gatk/guide/article?id=1204">here</a>.</li>
<li>It must be sorted in <strong>coordinate order</strong> (not by queryname and not &quot;unsorted&quot;).</li>
<li>It must list the <a href="http://www.broadinstitute.org/gatk/guide/article?id=6472">read groups</a> with sample names in the header.</li>
<li>Every read must belong to a read group.</li>
<li>The BAM file must pass Picard <a href="https://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile">ValidateSamFile</a> validation.</li>
</ul>
<p>See the <a href="http://samtools.github.io/hts-specs/">official BAM specification</a> for more information on what constitutes a valid BAM file.</p>
<hr />
<h3>4. What is the canonical ordering of human reference contigs in a BAM file?</h3>
<p>It depends on whether you're using the NCBI/GRC build 36/build 37 version of the human genome, or the UCSC hg18/hg19 version of the human genome.  While substantially equivalent, the naming conventions are different.  The canonical ordering of contigs for these genomes is as follows:</p>
<p>Human genome reference consortium standard ordering and names (b3x):
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT...</p>
<p>UCSC convention (hg1x):
chrM, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY...</p>
<hr />
<h3>5. How can I tell if my BAM file is sorted properly?</h3>
<p>The easiest way to do it is to download Samtools and run the following command to examine the header of your file:</p>
<pre><code class="pre_md">$ samtools view -H /path/to/my.bam
@HD     VN:1.0  GO:none SO:coordinate
@SQ     SN:1    LN:247249719
@SQ     SN:2    LN:242951149
@SQ     SN:3    LN:199501827
@SQ     SN:4    LN:191273063
@SQ     SN:5    LN:180857866
@SQ     SN:6    LN:170899992
@SQ     SN:7    LN:158821424
@SQ     SN:8    LN:146274826
@SQ     SN:9    LN:140273252
@SQ     SN:10   LN:135374737
@SQ     SN:11   LN:134452384
@SQ     SN:12   LN:132349534
@SQ     SN:13   LN:114142980
@SQ     SN:14   LN:106368585
@SQ     SN:15   LN:100338915
@SQ     SN:16   LN:88827254
@SQ     SN:17   LN:78774742
@SQ     SN:18   LN:76117153
@SQ     SN:19   LN:63811651
@SQ     SN:20   LN:62435964
@SQ     SN:21   LN:46944323
@SQ     SN:22   LN:49691432
@SQ     SN:X    LN:154913754
@SQ     SN:Y    LN:57772954
@SQ     SN:MT   LN:16571
@SQ     SN:NT_113887    LN:3994
...</code class="pre_md"></pre>
<p>If the order of the contigs here matches the contig ordering specified above, and the <code>SO:coordinate</code> flag appears in your header, then your contig and read ordering satisfies the GATK requirements.</p>
<hr />
<h3>6. My BAM file isn't sorted that way.  How can I fix it?</h3>
<p><a href="http://picard.sourceforge.net/">Picard</a> offers a tool called <a href="http://picard.sourceforge.net/command-line-overview.shtml#SortSam">SortSam</a>  that will sort a BAM file properly. A similar utility exists in Samtools, but we recommend the Picard tool because SortSam will also set a flag in the header that specifies that the file is correctly sorted, and this flag is necessary for the GATK to know it is safe to process the data.  Also, you can use the <a href="http://picard.sourceforge.net/command-line-overview.shtml">ReorderSam</a> command to make a BAM file SQ order match another reference sequence.</p>
<hr />
<h3>7. How can I tell if my BAM file has read group and sample information?</h3>
<p>A quick Unix command using Samtools will do the trick:</p>
<pre><code class="pre_md">$ samtools view -H /path/to/my.bam | grep '^@RG'
@RG ID:0    PL:solid    PU:Solid0044_20080829_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP   LB:Lib1 PI:2750 DT:2008-08-28T20:00:00-0400 SM:NA12414  CN:bcm
@RG ID:1    PL:solid    PU:0083_BCM_20080719_1_Pilot1_Ceph_12414_B_lib_1_2Kb_MP_Pilot1_Ceph_12414_B_lib_1_2Kb_MP    LB:Lib1 PI:2750 DT:2008-07-18T20:00:00-0400 SM:NA12414  CN:bcm
@RG ID:2    PL:LS454    PU:R_2008_10_02_06_06_12_FLX01080312_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
@RG ID:3    PL:LS454    PU:R_2008_10_02_06_07_08_rig19_retry    LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
@RG ID:4    PL:LS454    PU:R_2008_10_02_17_50_32_FLX03080339_retry  LB:HL#01_NA11881    PI:0    SM:NA11881  CN:454MSC
...</code class="pre_md"></pre>
<p>The presence of the <code>@RG</code> tags indicate the presence of read groups.  Each read group has a <code>SM</code> tag, indicating the sample from which the reads belonging to that read group originate.</p>
<p>In addition to the presence of a read group in the header, each read must belong to one and only one read group.  Given the following example reads,</p>
<pre><code class="pre_md">$ samtools view /path/to/my.bam | grep '^@RG'
EAS139_44:2:61:681:18781    35  1   1   0   51M =   9   59  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA B&lt;&gt;;==?=?&lt;==?=?=&gt;&gt;?&gt;&gt;&lt;=&lt;?=?8&lt;=?&gt;?&lt;:=?&gt;?&lt;==?=&gt;:;&lt;?:= RG:Z:4  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
EAS139_44:7:84:1300:7601    35  1   1   0   51M =   12  62  TAACCCTAAGCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA G&lt;&gt;;==?=?&amp;=&gt;?=?&lt;==?&gt;?&lt;&gt;&gt;?=?&lt;==?&gt;?&lt;==?&gt;?1==@&gt;?;&lt;=&gt;&lt;; RG:Z:3  MF:i:18 Aq:i:0  NM:i:1  UQ:i:5  H0:i:0  H1:i:85
EAS139_44:8:59:118:13881    35  1   1   0   51M =   2   52  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @&lt;&gt;;&lt;=?=?==&gt;?&gt;?&lt;==?=&gt;&lt;=&gt;?-?;=&gt;?:&gt;&lt;==?7?;&lt;&gt;?5?&lt;&lt;=&gt;:; RG:Z:1  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
EAS139_46:3:75:1326:2391    35  1   1   0   51M =   12  62  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA @&lt;&gt;==&gt;?&gt;@???B&gt;A&gt;?&gt;A?A&gt;??A?@&gt;?@A?@;??A&gt;@7&gt;?&gt;&gt;@:&gt;=@;@ RG:Z:0  MF:i:18 Aq:i:0  NM:i:0  UQ:i:0  H0:i:85 H1:i:31
...</code class="pre_md"></pre>
<p>membership in a read group is specified by the <code>RG:Z:*</code> tag.  For instance, the first read belongs to read group 4 (sample NA11881), while the last read shown here belongs to read group 0 (sample NA12414).</p>
<hr />
<h3>8. My BAM file doesn't have read group and sample information.  Do I really need it?</h3>
<p>Yes!  Many algorithms in the GATK need to know that certain reads were sequenced together on a specific lane, as they attempt to compensate for variability from one sequencing run to the next.  Others need to know that the data represents not just one, but many samples.  Without the read group and sample information, the GATK has no way of determining this critical information. You can use Picard's <a href="https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups">AddOrReplaceReadGroups</a> tool to add read group information.</p>
<hr />
<h3>11. What's the best way to create a subset of my BAM file containing only reads over a small interval?</h3>
<p>You can use the GATK to do the following:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -R reference.fasta -I full_input.bam -T PrintReads -L chr1:10-20 -o subset_input.bam</code class="pre_md"></pre>
<p>and you'll get a BAM file containing only reads overlapping those points.  This operation retains the complete BAM header from the full file (this was the reference aligned to, after all) so that the BAM remains easy to work with.  We routinely use these features for testing and high-performance analysis with the GATK.</p>