## (How to) Map reads to a reference with alternate contigs like GRCh38

http://gatkforums.broadinstitute.org/gatk/discussion/8017/how-to-map-reads-to-a-reference-with-alternate-contigs-like-grch38

<h4>Document is in <code>BETA</code>. It may be incomplete and/or inaccurate. Post suggestions to the <em>Comments</em> section and be sure to read about updates also within the <em>Comments</em> section.</h4>
<hr />
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/65/a0f09aad5f351a1322f7c1b19ec5d9.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/65/a0f09aad5f351a1322f7c1b19ec5d9.png" align="right" width="480" style="margin:5px 0px 5px 5px" /></a> This exploratory tutorial provides instructions and example data to map short reads to a reference genome with alternate haplotypes. Instructions are suitable for indexing and mapping reads to GRCh38. </p>
<p>► If you are unfamiliar with terms that describe reference genome components, or GRCh38 alternate haplotypes, take a few minutes to study the <em>Dictionary</em> entry <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">Reference Genome Components</a>.</p>
<p>► For an introduction to GRCh38, see <a href="https://software.broadinstitute.org/gatk/blog?id=8180">Blog#8180</a>.</p>
<p>Specifically, the tutorial uses BWA-MEM to index and map simulated reads for three samples to a mini-reference composed of a GRCh38 chromosome and alternate contig (<strong>sections 1–3</strong>). We align in an alternate contig aware (alt-aware) manner, which we also call alt-handling. This is the main focus of the tutorial. </p>
<p>The decision to align to a genome with alternate haplotypes has implications for variant calling. We discuss these in <strong>section 5</strong> using the callset generated with the optional tutorial steps outlined in <strong>section 4</strong>. Because we strategically placed a number of SNPs on the sequence used to simulate the reads, in both homologous and divergent regions, we can use the variant calls and their annotations to examine the implications of analysis approaches. To this end, the tutorial fast-forwards through pre-processing and calls variants for a trio of samples that represents the combinations of the two reference haplotypes (the PA and the ALT). This first workflow (<strong>tutorial_8017</strong>) is suitable for calling variants on the primary assembly but is insufficient for capturing variants on the alternate contigs.</p>
<p>For those who are interested in calling variants on the alternate contigs, we also present a second and a third workflow in <strong>section 6</strong>. The second workflow (<strong>tutorial_8017_toSE</strong>) takes the processed BAM from the first workflow, makes some adjustments to the reads to maximize their information, and calls variants on the alternate contig. This approach is suitable for calling on ~75% of the non-HLA alternate contigs or ~92% of loci with non-HLA alternate contigs (see <a href="https://us.v-cdn.net/5019796/uploads/FileUpload/fc/2834e6593da374296a205f33d117ac.png">table in section 6</a>). The third workflow (<strong>tutorial_8017_postalt</strong>) takes the alt-aware alignments from the first workflow and performs a postalt-processing step as well as the same adjustment from the second workflow. Postalt-processing uses the bwa-postalt.js javascript program that Heng Li provides as a companion to BWA. This allows for variant calling on all alternate contigs including HLA alternate contigs. </p>
<p>The tutorial ends by comparing the difference in call qualities from the multiple workflows <em>for the given example data</em> and discusses a few caveats of each approach.  </p>
<p><a name="top"></a></p>
<p>► The three workflows shown in the diagram above are available as <a href="https://software.broadinstitute.org/wdl/">WDL scripts</a> in our <a href="https://github.com/broadinstitute/wdl/tree/develop/scripts/tutorials/gatk">GATK Tutorials WDL scripts repository</a>.  </p>
<hr />
<h3>Jump to a section</h3>
<ol>
<li><a href="#1">Index the reference FASTA for use with BWA-MEM</a></li>
<li><a href="#2">Include the reference ALT index file</a>
☞ <a href="#2.1"><em>What happens if I forget the ALT index file?</em></a></li>
<li><a href="#3">Align reads with BWA-MEM</a>
☞ <a href="#3.1"><em>How can I tell if a BAM was aligned with alt-handling?</em></a>
☞ <a href="#3.2"><em>What is the <code>pa</code> tag?</em></a></li>
<li>(Optional) <a href="#4">Add read group information, preprocess to make a clean BAM and call variants</a></li>
<li><a href="#5">How can I tell  whether I should consider an alternate haplotype for a given sample?</a>
(5.1) <a href="#5.1">Discussion of variant calls for <strong>tutorial_8017</strong></a></li>
<li><a href="#6">My locus includes an alternate haplotype. How can I call variants on alt contigs?</a>
(6.1) <a href="#6.1">Variant calls for <strong>tutorial_8017_toSE</strong></a>
(6.2) <a href="#6.2">Variant calls for <strong>tutorial_8017_postalt</strong></a></li>
<li><a href="#7">Related resources</a></li>
</ol>
<h3>Tools involved</h3>
<ul>
<li>BWA v0.7.13 or later releases. The tutorial uses v0.7.15.
Download from <a href="https://sourceforge.net/projects/bio-bwa/files/">here</a> and see <a href="https://software.broadinstitute.org/gatk/documentation/article?id=2899">Tutorial#2899</a> for installation instructions.
The <code>bwa-postalt.js</code> script is within the <code>bwakit</code> folder.</li>
<li>Picard tools v2.5.0 or later releases. The tutorial uses v2.5.0.</li>
<li>Optional GATK tools. The tutorial uses v3.6.</li>
<li>Optional Samtools. The tutorial uses v1.3.1.</li>
<li>Optional <a href="https://www.gnu.org/software/gawk/">Gawk</a>, an <a href="https://en.wikipedia.org/wiki/AWK">AWK</a>-like tool that can interpret bitwise SAM flags. The tutorial uses v4.1.3.</li>
<li>Optional k8 Javascript shell. The tutorial uses v0.2.3 downloaded from <a href="https://github.com/attractivechaos/k8/releases/">here</a>.</li>
</ul>
<h3>Download example data</h3>
<p>Download tutorial_8017.tar.gz, either from the <a href="https://drive.google.com/open?id=0BzI1CyccGsZibnRtQjhaakxobEE">GoogleDrive</a> or from the <a href="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/tutorials/datasets">ftp site</a>. To access the ftp site, leave the password field blank. The data tarball contains the paired FASTQ reads files for three samples. It also contains a mini-reference <code>chr19_chr19_KI270866v1_alt.fasta</code> and corresponding <code>.dict</code> dictionary, <code>.fai</code> index and six BWA indices including the <code>.alt</code> index. The data tarball includes the output files from the workflow that we care most about. These are the aligned SAMs, processed and indexed BAMs and the final multisample VCF callsets from the three presented workflows.</p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/31/f1f2c77b6efbf9565700516b836914.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/31/f1f2c77b6efbf9565700516b836914.png" align="right" width="360" style="margin:10px 0px 10px 5px"/></a> The mini-reference contains two contigs subset from human GRCh38: <code>chr19</code> and <code>chr19_KI270866v1_alt</code>. The ALT contig corresponds to a diverged haplotype of chromosome 19. Specifically, it corresponds to chr19:34350807-34392977, which contains the <em>glucose-6-phosphate isomerase</em> or GPI gene. Part of the ALT contig introduces novel sequence that lacks a corresponding region in the primary assembly.</p>
<p>Using instructions in <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7859">Tutorial#7859</a>, we simulated paired 2x151 reads to derive three different sample reads that when aligned give roughly 35x coverage for the target primary locus. We derived the sequences from either the 43 kbp ALT contig (sample ALTALT), the corresponding 42 kbp region of the primary assembly (sample PAPA) or both (sample PAALT). Before simulating the reads, we introduced four SNPs to each contig sequence in a deliberate manner so that we can call variants.</p>
<p>► Alternatively, you may instead use the example input files and commands with the full <a href="https://software.broadinstitute.org/gatk/download/bundle">GRCh38 reference</a>. Results will be similar with a handful of reads mapping outside of the mini-reference regions.
<a name="1"></a></p>
<hr />
<h2>1. Index the reference FASTA for use with BWA-MEM</h2>
<p>Our example <code>chr19_chr19_KI270866v1_alt.fasta</code> reference already has <code>chr19_chr19_KI270866v1_alt.dict</code> dictionary and <code>chr19_chr19_KI270866v1_alt.fasta.fai</code> index files for use with Picard and GATK tools. BWA requires a different set of index files for alignment. The command below creates five of the six index files we need for alignment. The command calls the <code>index</code> function of BWA on the reference FASTA. </p>
<pre><code class="pre_md">bwa index chr19_chr19_KI270866v1_alt.fasta</code class="pre_md"></pre>
<p>This gives <code>.pac</code>, <code>.bwt</code>, <code>.ann</code>, <code>.amb</code> and <code>.sa</code> index files that all have the same <code>chr19_chr19_KI270866v1_alt.fasta</code> basename. Tools recognize index files within the same directory by their identical basename. In the case of BWA, it uses the basename preceding the <code>.fasta</code> suffix and searches for the index file, e.g. with <code>.bwt</code> suffix or <code>.64.bwt</code> suffix. Depending on which of the two choices it finds, it looks for the same suffix for the other index files, e.g. <code>.alt</code> or <code>.64.alt</code>. Lack of a matching <code>.alt</code> index file will cause BWA to map reads without alt-handling. More on this next.  </p>
<p>Note that the <code>.64.</code> part is an explicit indication that index files were generated with version 0.6 or later of BWA and are the 64-bit indices (as opposed to files generated by earlier versions, which were 32-bit). This <code>.64.</code> signifier can be added automatically by adding <code>-6</code> to the <code>bwa index</code> command. </p>
<p><a name="2"></a>
<a href="#top">back to top</a></p>
<hr />
<h2>2. Include the reference ALT index file</h2>
<p>Be sure to place the tutorial's mini-ALT index file <code>chr19_chr19_KI270866v1_alt.fasta.alt</code> with the other index files. Also, if it does not already match, change the file basename to match. This is the sixth index file we need for alignment. BWA-MEM uses this file to prioritize primary assembly alignments for reads that can map to both the primary assembly and an alternate contig. See <a href="https://github.com/lh3/bwa/blob/master/README-alt.md">BWA documentation</a> for details.</p>
<ul>
<li>As of this writing (August 8, 2016), the SAM format ALT index file for GRCh38 is available only in the <a href="https://sourceforge.net/projects/bio-bwa/files/bwakit/">x86_64-linux bwakit download</a> as stated in this <a href="https://github.com/lh3/bwa/tree/master/bwakit">bwakit README</a>. The <code>hs38DH.fa.alt</code> file is in the <code>resource-GRCh38</code> folder.</li>
<li>In addition to <em>mapped</em> alternate contig records, the ALT index also contains decoy contig records as <em>unmapped</em> SAM records. This is relevant to the postalt-processing we discuss in <a href="#6.2">section 6.2</a>. As such, the postalt-processing in <strong>section 6</strong> also requires the ALT index. </li>
</ul>
<p>For the tutorial, we subset from <code>hs38DH.fa.alt</code> to create a mini-ALT index, <code>chr19_chr19_KI270866v1_alt.fasta.alt</code>. Its contents are shown below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/0b/8b7f95e734b3d9faa45ba8f8deea28.png" />
<p><a name="2.1"></a></p>
<p>The record aligns the <code>chr19_KI270866v1_alt</code> contig to the <code>chr19</code> locus starting at position 34,350,807 and uses CIGAR string nomenclature to indicate the pairwise structure. To interpret the CIGAR string, think of the primary assembly as the reference and the ALT contig sequence as the read. For example, the <code>11307M</code> at the start indicates 11,307 corresponding sequence bases, either matches or mismatches. The <code>935S</code> at the end indicates a 935 base softclip for the ALT contig sequence that lacks corresponding sequence in the primary assembly. This is a region that we consider highly divergent or novel. Finally, notice the <code>NM</code> tag that notes the edit distance to the reference.</p>
<h3>☞ What happens if I forget the ALT index file?</h3>
<p>If you omit the ALT index file from the reference, or if its naming structure mismatches the other indexes, then your alignments will be equivalent to the results you would obtain if you run BWA-MEM with the <code>-j</code> option. The next section gives an example of what this looks like.</p>
<p><a name="3"></a>
<a href="#top">back to top</a></p>
<hr />
<h2>3. Align reads with BWA-MEM</h2>
<p>The command below uses an alt-aware version of BWA and maps reads using BWA's <em>maximal exact match</em> (MEM) option. Because the ALT index file is present, the tool prioritizes mapping to the primary assembly over ALT contigs. In the command, the tutorial's <code>chr19_chr19_KI270866v1_alt.fasta</code> serves as reference; one FASTQ holds the forward reads and the other holds the reverse reads.  </p>
<pre><code class="pre_md">bwa mem chr19_chr19_KI270866v1_alt.fasta 8017_read1.fq 8017_read2.fq &gt; 8017_bwamem.sam</code class="pre_md"></pre>
<p>The resulting file <code>8017_bwamem.sam</code> contains aligned read records. </p>
<p><a name="3.1"></a></p>
<ul>
<li>BWA preferentially maps to the primary assembly any reads that can align equally well to the primary assembly or the ALT contigs as well as any reads that it can reasonably align to the primary assembly even if it aligns better to an ALT contig. Preference is given by the <em>primary</em> alignment record status, i.e. not <em>secondary</em> and not <em>supplementary</em>. BWA takes the reads that it cannot map to the primary assembly and attempts to map them to the alternate contigs. If a read can map to an alternate contig, then it is mapped to the alternate contig as a <em>primary</em> alignment. For those reads that can map to both and align better to the ALT contig, the tool flags the ALT contig alignment record as <em>supplementary</em> (0x800). This is what we call alt-aware mapping or alt-handling.</li>
<li>Adding the <code>-j</code> option to the command disables the alt-handling. Reads that can map multiply are given low or zero MAPQ scores.</li>
</ul>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/38/98aa1f4e0468b7fc8106a6bcc600c5.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/38/98aa1f4e0468b7fc8106a6bcc600c5.png" align="right" width="450" style="margin:20px 0px 5px 5px"/></a> </p>
<h3>☞ How can I tell if a BAM was aligned with alt-handling?</h3>
<p>There are two approaches to this question.</p>
<p>First, you can view the alignments on IGV and compare primary assembly loci with their alternate contigs. The IGV screenshots to the right show how BWA maps reads with (top) or without (bottom) alt-handling. </p>
<p>Second, you can check the alignment SAM. Of two tags that indicate alt-aware alignment, one will persist after preprocessing only if the sample has reads that can map to alternate contigs. The first tag, the <code>AH</code> tag, is in the BAM header section of the alignment file, and is absent after any merging step, e.g. merging with MergeBamAlignment. The second tag, the <code>pa</code> tag, is present for reads that the aligner alt-handles. If a sample does not contain any reads that map equally or preferentially to alternate contigs, then this tag may be absent in a BAM even if the alignments were mapped in an alt-aware manner.</p>
<p>Here are three headers for comparison where only one <em>indicates</em> alt-aware alignment.</p>
<p><strong>File header for alt-aware alignment. We use this type of alignment in the tutorial.</strong>
Each alternate contig's <code>@SQ</code> line in the header will have an <code>AH:*</code> tag to indicate alternate contig handling for that contig. This marking is based on the alternate contig being listed in the <code>.alt</code> index file and alt-aware alignment. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/68/286e35c8b4b91648ab2a632d41d60d.png" />
<p><strong>File header for <code>-j</code> alignment (alt-handling disabled) for example purposes. We do not perform this type of alignment in the tutorial.</strong>
Notice the absence of any special tags in the header.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/cf/e945f8f24e7b0ecbf37596fb05063b.png" />
<p><a name="3.2"></a></p>
<p><strong>File header for alt-aware alignment after merging with MergeBamAlignment. We use this step in the next section.</strong>
Again, notice the absence of any special tags in the header.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/5d/b06974ff127d9fa74ac71071d87b7f.png" />
<h3>☞ What is the <code>pa</code> tag?</h3>
<p>For BWA v0.7.15, but not v0.7.13, ALT loci alignment records that can align to both the primary assembly and alternate contig(s) will have a <code>pa</code> tag on the primary assembly alignment. For example, read <code>chr19_KI270866v1_alt_4hetvars_26518_27047_0:0:0_0:0:0_931</code> of the ALTALT sample has five alignment records only three of which have the <code>pa</code> tag as shown below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/f7/a79b13a178f886056a793556ccbf85.png" />
<p>A brief description of each of the five alignments, in order:</p>
<ol>
<li>First in pair, primary alignment on the <em>primary assembly</em>; AS=146, <strong>pa=0.967</strong> </li>
<li>First in pair, supplementary alignment on the <em>alternate contig</em>; AS=151</li>
<li>Second in pair, primary alignment on the <em>primary assembly</em>; AS=120; <strong>pa=0.795</strong> </li>
<li>Second in pair, supplementary alignment on the <em>primary assembly</em>; AS=54; <strong>pa=0.358</strong> </li>
<li>Second in pair, supplementary alignment on the <em>alternate contig</em>; AS=151 </li>
</ol>
<p>The <code>pa</code> tag measures how much better a read aligns to its best alternate contig alignment versus its primary assembly (pa) alignment. Specifically, it is the ratio of the primary assembly alignment score over the highest alternate contig alignment score. In our example we have primary assembly alignment scores of 146, 120 and 54 and alternate contig alignment scores of 151 and again 151. This gives us three different <code>pa</code> scores that tag the primary assembly alignments: 146/151=0.967, 120/151=0.795 and 54/151=0.358. </p>
<p>In our tutorial's workflow, MergeBamAlignment may either change an alignment's <code>pa</code> score or add a previously unassigned <code>pa</code> score to an alignment. The result of this is summarized as follows for the same alignments. </p>
<ol>
<li>pa=0.967 --MergeBamAlignment--&gt; same</li>
<li>none --MergeBamAlignment--&gt; assigns pa=0.967</li>
<li>pa=0.795 --MergeBamAlignment--&gt; same</li>
<li>pa=0.358 --MergeBamAlignment--&gt; <strong>replaces with pa=0.795</strong></li>
<li>none --MergeBamAlignment--&gt; assigns pa=0.795</li>
</ol>
<p>If you want to retain the BWA-assigned <code>pa</code> scores, then add the following options to the workflow commands in <strong>section 4</strong>.</p>
<ul>
<li>For RevertSam, add <code>ATTRIBUTE_TO_CLEAR=pa</code>.</li>
<li>For MergeBamAlignment, add <code>ATTRIBUTES_TO_RETAIN=pa</code>.</li>
</ul>
<p>In our sample set, after BWA-MEM alignment ALTALT has 1412 <code>pa</code>-tagged alignment records, PAALT has 805 <code>pa</code>-tagged alignment records and PAPA has zero <code>pa</code>-tagged records. </p>
<p><a name="4"></a>
<a href="#top">back to top</a></p>
<hr />
<h2>4. Add read group information, preprocess to make a clean BAM and call variants</h2>
<p>The initial alignment file is missing read group information. One way to add that information, which we use in production, is to use <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeBamAlignment">MergeBamAlignment</a>. MergeBamAlignment adds back read group information contained in an unaligned BAM and adjusts meta information to produce a clean BAM ready for pre-processing (see <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6483">Tutorial#6483</a> for details on our use of MergeBamAlignment). Given the focus here is to showcase BWA-MEM's alt-handling, we refrain from going into the details of all this additional processing. They follow, with some variation, the PairedEndSingleSampleWf pipeline detailed <a href="https://github.com/broadinstitute/wdl/blob/develop/scripts/broad_pipelines/PairedSingleSampleWf_160720.md">here</a>.</p>
<p>Remember these are simulated reads with simulated base qualities. We simulated the reads in a manner that only introduces the planned mismatches, without any errors. Coverage is good at roughly 35x. All of the base qualities for all of the reads are at <code>I</code>, which is, according to <a href="https://en.wikipedia.org/wiki/FASTQ_format">this page</a> and <a href="http://broadinstitute.github.io/picard/explain-qualities.html">this site</a>, an excellent base quality score equivalent to a Sanger Phred+33 score of 40. We can therefore skip base quality score recalibration (BQSR) since the reads are simulated and the dataset is not large enough for recalibration anyway. </p>
<p>Here are the commands to obtain a final multisample variant callset. The commands are given for one of the samples. Process each of the three samples independently in the same manner [4.1–4.6] until the last GenotypeGVCFs command [4.7].</p>
<p>[4.1] Create unmapped uBAM</p>
<pre><code>java -jar picard.jar RevertSam \
    I=altalt_bwamem.sam O=altalt_u.bam \
    ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA</code></pre>
<p>[4.2] Add read group information to uBAM</p>
<pre><code>java -jar picard.jar AddOrReplaceReadGroups \
    I=altalt_u.bam O=altalt_rg.bam \
    RGID=altalt RGSM=altalt RGLB=wgsim RGPU=shlee RGPL=illumina</code></pre>
<p>[4.3] Merge uBAM with aligned BAM</p>
<pre><code>java -jar picard.jar MergeBamAlignment \
    ALIGNED=altalt_bwamem.sam UNMAPPED=altalt_rg.bam O=altalt_m.bam \
    R=chr19_chr19_KI270866v1_alt.fasta \
    SORT_ORDER=unsorted CLIP_ADAPTERS=false \
    ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    UNMAP_CONTAMINANT_READS=false \
    ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA</code></pre>
<p>[4.4] Flag duplicate reads</p>
<pre><code>java -jar picard.jar MarkDuplicates \
    INPUT=altalt_m.bam OUTPUT=altalt_md.bam METRICS_FILE=altalt_md.bam.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname </code></pre>
<p>[4.5] Coordinate sort, fix <code>NM</code> and <code>UQ</code> tags and index for clean BAM
As of Picard v2.7.0, released October 17, 2016, <strong>SetNmAndUqTags</strong> is no longer available. Use <strong>SetNmMdAndUqTags</strong> instead.</p>
<pre><code>set -o pipefail
java -jar picard.jar SortSam \
    INPUT=altalt_md.bam OUTPUT=/dev/stdout SORT_ORDER=coordinate | \
    java -jar $PICARD SetNmAndUqTags \
    INPUT=/dev/stdin OUTPUT=altalt_snaut.bam \
    CREATE_INDEX=true R=chr19_chr19_KI270866v1_alt.fasta</code></pre>
<p>[4.6] Call SNP and indel variants in <em>emit reference confidence</em> (ERC) mode per sample using HaplotypeCaller</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R chr19_chr19_KI270866v1_alt.fasta \
    -o altalt.g.vcf -I altalt_snaut.bam \
    -ERC GVCF --max_alternate_alleles 3 --read_filter OverclippedRead \
    --emitDroppedReads -bamout altalt_hc.bam</code></pre>
<p>[4.7] Call genotypes on three samples</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -R chr19_chr19_KI270866v1_alt.fasta -o multisample.vcf \
    --variant altalt.g.vcf --variant altpa.g.vcf --variant papa.g.vcf </code></pre>
<p>The <code>altalt_snaut.bam</code>, HaplotypeCaller's <code>altalt_hc.bam</code> and the multisample <code>multisample.vcf</code> are ready for viewing on IGV.</p>
<p>Before getting into the results in the next section, we have minor comments on two filtering options. </p>
<p>In our tutorial workflows, we turn off MergeBamAlignment's <code>UNMAP_CONTAMINANT_READS</code> option. If set to true, 68 reads become unmapped for PAPA and 40 reads become unmapped for PAALT. These unmapped reads are those reads caught by the <code>UNMAP_CONTAMINANT_READS</code> filter <em>and their mates</em>. MergeBamAlignment defines contaminant reads as those alignments that are overclipped, i.e. that are softclipped on both ends, and that align with less than 32 bases. Changing the <code>MIN_UNCLIPPED_BASES</code> option from the default of 32 to 22 and 23 restores all of these reads for PAPA and PAALT, respectively. Contaminants are obviously absent for these simulated reads. And so we set <code>UNMAP_CONTAMINANT_READS</code> to false to disable this filtering.</p>
<p>HaplotypeCaller's <code>--read_filter OverclippedRead</code> option similarly looks for both-end-softclipped alignments, then filters reads aligning with less than 30 bases. The difference is that HaplotypeCaller only excludes the overclipped alignments from its calling and does not remove mapping information nor does it act on the mate of the filtered alignment. Thus, we keep this read filter for the first workflow. However, for the second and third workflows in <a href="#6">section 6</a>, <strong>tutorial_8017_toSE</strong> and <strong>tutorial_8017_postalt</strong>, we omit the <code>--read_filter Overclipped</code> option from the HaplotypeCaller command. We also omit the <code>--max_alternate_alleles 3</code> option for simplicity.</p>
<p><a name="5"></a>
<a href="#top">back to top</a></p>
<hr />
<h2>5. How can I tell  whether I should consider an alternate haplotype?</h2>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/30/3d6bbf9a4c67674e1ebea0308bdd3f.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/30/3d6bbf9a4c67674e1ebea0308bdd3f.png" align="right" width="450" style="margin:20px 0px 5px 5px"/></a> We consider this question only for our GPI locus, a locus we know has an alternate contig in the reference. Here we use the term <em>locus</em> in its biological sense to refer to a contiguous genomic region of interest. The three samples give the alignment and coverage profiles shown on the right.  </p>
<p>What is immediately apparent from the IGV screenshot is that the scenarios that include the alternate haplotype give a distinct pattern of variant sites to the primary assembly much like a fingerprint. These variants are predominantly heterozygous or homozygous. Looking closely at the 3' region of the locus, we see some alignment coverage anomalies that also show a distinct pattern. The coverage in some of the highly diverged region in the primary assembly drops while in others it increases. If we look at the origin of simulated reads in one of the excess coverage regions, we see that they are from two different regions of the alternate contig that suggests duplicated sequence segments within the alternate locus.</p>
<p>The variation pattern and coverage anomalies on the primary locus suggest an alternate haplotype may be present for the locus. We can then confirm the presence of aligned reads, both supplementary and primary, on the alternate locus. Furthermore, if we count the alignment records for each region, e.g. using <code>samtools idxstats</code>, we see the following metrics. </p>
<pre><code>                        ALT/ALT     PA/ALT     PA/PA   
chr19                     10005      10006     10000     
chr19_KI270866v1_alt       1407        799         0      </code></pre>
<p><a name="5.1"></a></p>
<p>The number of alignments on the alternate locus increases proportionately with alternate contig dosage. All of these factors together suggest that the sample presents an alternate haplotype. </p>
<h3>5.1 Discussion of variant calls for tutorial_8017</h3>
<p>The three-sample variant callset gives 54 sites on the primary locus and two additional on the alternate locus for 56 variant sites. All of the eight SNP alleles we introduced are called, with six called on the primary assembly and two called on the alternate contig. Of the 15 expected genotype calls, four are incorrect. Namely, four PAALT calls that ought to be heterozygous are called homozygous variant. These are two each on the primary assembly and on the alternate contig in the region that is highly divergent.</p>
<p>► Our production pipelines use genomic intervals lists that exclude GRCh38 alternate contigs from <em>variant calling</em>. That is, variant calling is performed only for contigs of the primary assembly. This calling on even just the primary assembly of GRCh38 brings improvements to analysis results over previous assemblies. For example, if we align and call variants for our simulated reads on GRCh37, we call 50 variant sites with identical QUAL scores to the equivalent calls in our GRCh38 callset. However, this GRCh37 callset is missing six variant calls compared to the GRCh38 callset for the 42 kb locus: the two variant sites on the alternate contig and <em>four variant sites on the primary assembly</em>. </p>
<p>Consider the example variants on the primary locus. The variant calls from the primary assembly include 32 variant sites that are strictly homozygous variant in ALTALT <em>and</em> heterozygous variant in PAALT. The callset represents only those reads from the ALT <em>that can be mapped to the primary assembly</em>.</p>
<p>In contrast, the two variants in regions whose reads <em>can only map to the alternate contig</em> are absent from the primary assembly callset. For this simulated dataset, the primary alignments present on the alternate contig provide enough supporting reads that allow HaplotypeCaller to call the two variants. However, these variant calls have <em>lower-quality annotation metrics</em> than for those simulated in an equal manner on the primary assembly. We will get into why this is in <strong>section 6</strong>. </p>
<p>Additionally, for our PAALT sample that is heterozygous for an alternate haplotype, the genotype calls in the highly divergent regions are inaccurate. These are called homozygous variant on the primary assembly and on the alternate contig when in fact they are heterozygous variant. These calls have lower genotype scores <code>GQ</code> as well as lower allele depth <code>AD</code> and coverage <code>DP</code>. The table below shows the variant calls for the introduced SNP sites. In blue are the genotype calls that should be heterozygous variant but are instead called homozygous variant.
<a href="https://us.v-cdn.net/5019796/uploads/FileUpload/0a/bc1a0d986fcf6087058c9ce46551bf.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/0a/bc1a0d986fcf6087058c9ce46551bf.png" align="" width="" style="margin:10px 0px 5px 0px"/></a></p>
<p>Here is a command to select out the intentional variant sites that uses <a href="https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php">SelectVariants</a>:</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T SelectVariants \
    -R chr19_chr19_KI270866v1_alt.fasta \
    -V multisample.vcf -o multisample_selectvariants.vcf \
    -L chr19:34,383,500 -L chr19:34,389,485 -L chr19:34,391,800 -L chr19:34,392,600 \
    -L chr19_KI270866v1_alt:32,700 -L chr19_KI270866v1_alt:38,700 \
    -L chr19_KI270866v1_alt:41,700 -L chr19_KI270866v1_alt:42,700 \
    -L chr19:34,383,486 -L chr19_KI270866v1_alt:32,714 </code></pre>
<p><a name="6"></a>
<a href="#top">back to top</a></p>
<hr />
<h2>6. My locus includes an alternate haplotype. How can I call variants on alt contigs?</h2>
<p>If you want to call variants on alternate contigs, consider additional data processing that overcome the following problems.</p>
<ul>
<li>Loss of alignments from filtering of overclipped reads.</li>
<li>HaplotypeCaller's filtering of alignments whose mates map to another contig. Alt-handling produces many of these types of reads on the alternate contigs.</li>
<li>Zero MAPQ scores for alignments that map to two or more alternate contigs. HaplotypeCaller excludes these types of reads from contributing to evidence for variation. </li>
</ul>
<p>Let us talk about these in more detail.</p>
<p>Ideally, if we are interested in alternate haplotypes, then we would have ensured we were using the most up-to-date analysis reference genome sequence with the latest patch fixes. Also, whatever approach we take to align and preprocess alignments, if we filter any reads as putative contaminants, e.g. with MergeBamAlignment's option to unmap cross-species contamination, then at this point we would want to fish back into the unmapped reads pool and pull out those reads. Specifically, these would have an <code>SA</code> tag indicating mapping to the alternate contig of interest and an <code>FT</code> tag indicating the reason for unmapping was because MergeBamAlignment's <code>UNMAP_CONTAMINANT_READS</code> option identified them as cross-species contamination. Similarly, we want to make sure not to include HaplotypeCaller's <code>--read_filter OverclippedRead</code> option that we use in the first workflow. </p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/28/f8c0855ffef62382f0e96e53a82977.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/28/f8c0855ffef62382f0e96e53a82977.png" align="right" width="480" style="margin:5px 0px 5px 5px" /></a> As <strong>section 5.1</strong> shows, variant calls on the alternate contig are of low quality--they have roughly an order of magnitude lower QUAL scores than what should be equivalent variant calls on the primary assembly. </p>
<p>For this exploratory tutorial, we are interested in calling the introduced SNPs with equivalent annotation metrics. Whether they are called on the primary assembly or the alternate contig and whether they are called homozygous variant or heterozygous--let's say these are less important, especially given pinning certain variants from highly homologous regions to one of the loci is nigh impossible with our short reads. To this end, we will use the second workflow shown in the workflows diagram. However, because this solution is limited, we present a third workflow as well.  </p>
<p>► We present these workflows solely for exploratory purposes. They do not represent any production workflows.</p>
<p><strong>Tutorial_8017_toSE</strong> uses the processed BAM from our first workflow and allows for calling on singular alternate contigs. That is, the workflow is suitable for calling on alternate contigs of loci with only a single alternate contig like our GPI locus. <strong>Tutorial_8017_postalt</strong> uses the aligned SAM from the first workflow before processing, and requires separate processing before calling. This third workflow allows for calling on all alternate contigs, even on HLA loci that have numerous contigs per primary locus. However, the callset will not be parsimonious. That is, each alternate contig will greedily represent alignments and it is possible the same variant is called for all the alternate loci for a given primary locus as well as on the primary locus. It is up to the analyst to figure out what to do with the resulting calls. </p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/fc/2834e6593da374296a205f33d117ac.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/fc/2834e6593da374296a205f33d117ac.png" align="left" width="450" style="margin: 0px 5px 5px 0px"/></a> The reason for the divide in these two workflows is in the way BWA assigns mapping quality scores (MAPQ) to multimapping reads. Postalt-processing becomes necessary for loci with two or more alternate contigs because the shared alignments between the primary locus and alternate loci will have zero MAPQ scores. Postalt-processing gives non-zero MAPQ scores to the alignment records. The table presents the frequencies of GRCh38 non-HLA alternate contigs per primary locus. It appears that ~75% of non-HLA alternate contigs are singular to ~92% of primary loci with non-HLA alternate contigs. In terms of bases on the primary assembly, of the ~75 megabases that have alternate contigs, ~64 megabases (85%) have singular non-HLA alternate contigs and ~11 megabases (15%) have multiple non-HLA alternate contigs per locus. Our tutorial's example locus falls under this majority.</p>
<p><a name="6.1"></a></p>
<p>In both alt-aware mapping and postalt-processing, alternate contig alignments have a predominance of mates that map back to the primary assembly. HaplotypeCaller, for good reason, filters reads whose mates map to a different contig. However, we know that GRCh38 <em>artificially</em> represents alternate haplotypes as separate contigs and BWA-MEM <em>intentionally</em> maps these mates back to the primary locus. For comparable calls on alternate contigs, we need to include these alignments in calling. To this end, we have devised a temporary workaround. </p>
<h3>6.1 Variant calls for <strong>tutorial_8017_toSE</strong></h3>
<p>Here we are only aiming for <em>equivalent calls</em> with similar annotation values for the two variants that are called on the alternate contig. For the solution that we will outline, here are the results.  </p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/81/62cd9b4e7c710dc94ea6c4bdb45db9.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/81/62cd9b4e7c710dc94ea6c4bdb45db9.png" align="" width="" style="margin:10px 0px 5px 0px"/></a></p>
<p>Including the mate-mapped-to-other-contig alignments bolsters the variant call qualities for the two SNPs HaplotypeCaller calls on the alternate locus. We see the <code>AD</code> allele depths much improved for ALTALT and PAALT. Corresponding to the increase in reads, the <code>GQ</code> genotype quality and the QUAL score (highlighted in red) indicate higher qualities. For example, the QUAL scores increase from 332 and 289 to 2166 and 1764, respectively. We also see that one of the genotype calls changes. For sample ALTALT, we see a previous <em>no call</em> is now a homozygous reference call (highlighted in blue). This hom-ref call is further from the truth than not having a call as the ALTALT sample should not have coverage for this region in the primary assembly.</p>
<p>For our example data, <strong>tutorial_8017</strong>'s callset subset for the primary assembly and <strong>tutorial_8017_toSE</strong>'s callset subset for the alternate contigs together appear to make for a better callset. </p>
<p>What solution did we apply? As the workflow's name <em>toSE</em> implies, this approach converts paired reads to single end reads. Specifically, this approach takes the processed and coordinate-sorted BAM from the first workflow and removes the 0x1 <em>paired</em> flag from the alignments. Removing the 0x1 flag from the reads allows HaplotypeCaller to consider alignments whose mates map to a different contig. We accomplish this using a modified script of that presented in <em>Biostars</em> post <a href="https://www.biostars.org/p/106668/"><a href="https://www.biostars.org/p/106668/">https://www.biostars.org/p/106668/</a></a>, indexing with Samtools and then calling with HaplotypeCaller as follows. Note this workaround creates an invalid BAM according to ValidateSamFile. Also, another caveat is that because HaplotypeCaller uses softclipped sequences, any overlapping regions of read pairs will count twice towards variation instead of once. Thus, this step may lead to overconfident calls in such regions.  </p>
<p>Remove the 0x1 bitwise flag from alignments</p>
<pre><code>samtools view -h altalt_snaut.bam | gawk '{printf "%s\t", $1; if(and($2,0x1))
{t=$2-0x1}else{t=$2}; printf "%s\t" , t; for (i=3; i&lt;NF; i++){printf "%s\t", $i} ; 
printf "%s\n",$NF}'| samtools view -Sb - &gt; altalt_se.bam</code></pre>
<p>Index the resulting BAM</p>
<pre><code>samtools index altalt_se.bam</code></pre>
<p>Call variants in <code>-ERC GVCF</code> mode with HaplotypeCaller for each sample</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R chr19_chr19_KI270866v1_alt.fasta \
    -I altalt_se.bam -o altalt_hc.g.vcf \
    -ERC GVCF --emitDroppedReads -bamout altalt_hc.bam</code></pre>
<p><a name="6.2"></a></p>
<p>Finally, use GenotypeGVCFs as shown in <a href="#4">section 4</a>'s command [4.7] for a multisample variant callset. Tutorial_8017_toSE calls 68 variant sites--66 on the primary assembly and two on the alternate contig.</p>
<h3>6.2 Variant calls for <strong>tutorial_8017_postalt</strong></h3>
<p>BWA's postalt-processing requires the query-grouped output of BWA-MEM. Piping an alignment step with postalt-processing is possible. However, to be able to compare variant calls from an identical alignment, we present the postalt-processing as an <em>add-on</em> workflow that takes the alignment from the first workflow. </p>
<p>The command uses the <code>bwa-postalt.js</code> script, which we run through <code>k8</code>, a Javascript execution shell. It then lists the ALT index, the aligned SAM <code>altalt.sam</code> and names the resulting file <code>&gt; altalt_postalt.sam</code>. </p>
<pre><code>k8 bwa-postalt.js \
    chr19_chr19_KI270866v1_alt.fasta.alt \
    altalt.sam &gt; altalt_postalt.sam</code></pre>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/a5/3522324635aec94071a9ff688a4aa6.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a5/3522324635aec94071a9ff688a4aa6.png" align="right" width="450" style="margin:10px 0px 5px 5px"/></a> The resulting postalt-processed SAM, <code>altalt_postalt.sam</code>, undergoes the same processing as the first workflow (commands 4.1 through 4.7) except that (i) we omit <code>--max_alternate_alleles 3</code> and <code>--read_filter OverclippedRead</code> options for the HaplotypeCaller command like we did in <strong>section 6.1</strong> and (ii) we perform the 0x1 flag removal step from <strong>section 6.1</strong>.</p>
<p>The effect of this postalt-processing is immediately apparent in the IGV screenshots. Previously empty regions are now filled with alignments. Look closely in the highly divergent region of the primary locus. Do you notice a change, albeit subtle, before and after postalt-processing for samples ALTALT and PAALT?</p>
<p>These alignments give the calls below for our SNP sites of interest. Here, notice calls are made for more sites--on the equivalent site if present in addition to the design site (highlighted in the first two columns). For the three pairs of sites that can be called on either the primary locus or alternate contig, the variant site QUALs, the INFO field annotation metrics and the sample level annotation values are identical for each pair.</p>
<p><a href="https://us.v-cdn.net/5019796/uploads/FileUpload/c4/5f07b27798374175ba40f970e77a62.png"><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c4/5f07b27798374175ba40f970e77a62.png" align="" width="" style="margin:10px 0px 5px 0px"/></a> </p>
<p>Postalt-processing lowers the MAPQ of primary locus alignments in the highly divergent region that map better to the alt locus. You can see this as a subtle change in the IGV screenshot. After postalt-processing we see an increase in white zero MAPQ reads in the highly divergent region of the primary locus for ALTALT and PAALT. For ALTALT, this effectively cleans up the variant calls in this region at chr19:34,391,800 and chr19:34,392,600. Previously for ALTALT, these calls contained some reads: 4 and 25 for the first workflow and 0 and 28 for the second workflow. After postalt-processing, no reads are considered in this region giving us <code>./.:0,0:0:.:0,0,0</code> calls for both sites. </p>
<p>What we omit from examination are the effects of postalt-processing on decoy contig alignments. Namely, if an alignment on the primary assembly aligns better on a decoy contig, then postalt-processing discounts the alignment on the primary assembly by assigning it a zero MAPQ score.</p>
<p>To wrap up, here are the number of variant sites called for the three workflows. As you can see, this last workflow calls the most variants at 95 variant sites, with 62 on the primary assembly and 33 on the alternate contig.</p>
<pre><code>Workflow                total    on primary assembly    on alternate contig
tutorial_8017           56       54                      2
tutorial_8017_toSE      68       66                      2
tutorial_8017_postalt   95       62                     33</code></pre>
<p><a name="7"></a>
<a href="#top">back to top</a></p>
<hr />
<h3>7. Related resources</h3>
<ul>
<li>For WDL scripts of the workflows represented in this tutorial, see the <a href="https://github.com/broadinstitute/wdl/tree/develop/scripts/tutorials/gatk">GATK WDL scripts repository</a>. </li>
<li>To revert an aligned BAM to unaligned BAM, see <strong>Section B</strong> of <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6484">Tutorial#6484</a>.</li>
<li>To simulate reads from a reference contig, see <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7859">Tutorial#7859</a>.</li>
<li><em>Dictionary</em> entry <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">Reference Genome Components</a> reviews terminology that describe reference genome components.</li>
<li>The <a href="https://software.broadinstitute.org/gatk/download/bundle">GATK resource bundle</a> provides an analysis set GRCh38 reference FASTA as well as several other related resource files. </li>
<li>As of this writing (August 8, 2016), the SAM format ALT index file for GRCh38 is available only in the <a href="https://sourceforge.net/projects/bio-bwa/files/bwakit/">x86_64-linux bwakit download</a> as stated in this <a href="https://github.com/lh3/bwa/tree/master/bwakit">bwakit README</a>. The <code>hs38DH.fa.alt</code> file is in the <code>resource-GRCh38</code> folder. Rename this file's basename to match that of the corresponding reference FASTA.</li>
<li>For more details on MergeBamAlignment features, see <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6483#step3C">Section 3C</a> of <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6483">Tutorial#6483</a>.</li>
<li>For details on the PairedEndSingleSampleWorkflow that uses GRCh38, see <a href="https://github.com/broadinstitute/wdl/blob/develop/scripts/broad_pipelines/PairedSingleSampleWf_160720.md">here</a>.</li>
<li>See <a href="https://samtools.github.io/hts-specs">here</a> for VCF specifications.  </li>
</ul>
<p><a href="#top">back to top</a></p>
<hr />