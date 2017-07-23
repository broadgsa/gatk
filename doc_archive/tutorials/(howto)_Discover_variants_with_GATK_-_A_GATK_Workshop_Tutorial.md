## (howto) Discover variants with GATK - A GATK Workshop Tutorial

http://gatkforums.broadinstitute.org/gatk/discussion/7869/howto-discover-variants-with-gatk-a-gatk-workshop-tutorial

<h2>GATK TUTORIAL :: Variant Discovery :: Worksheet</h2>
<p><strong>June 2016 - GATK 3.6</strong></p>
<p>This tutorial covers material taught at GATK workshops, and focuses on key steps of the GATK Best Practices for Germline SNP and Indel Discovery in Whole Genomes and Exomes. If you aren't already, please set up your computer using the <a href="https://www.broadinstitute.org/gatk/guide/article?id=7098">workshop-specific installation instructions</a>. You can find additional background information relevant to this tutorial in the <a href="https://www.broadinstitute.org/gatk/guide/article?id=7870">Variant Discovery Appendix</a>.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/86/c38782834bee194dfc19ce3c37f40a.png" align=left width="200" hspace="10" vspace="10" />
<p>Our main purpose is to <strong>demonstrate an effective workflow for calling germline SNPs and indels</strong> in cohorts of multiple samples. This workflow can be applied to <strong>whole genomes</strong> as well as <strong>exomes</strong> and other targeted sequencing datasets. </p>
<p>We’ll start by examining the <strong>differences between data types</strong> (whole genomes, exomes and RNAseq) to highlight the properties of the data that influence what we need to do to analyze it as well as what we can expect to get out of it. </p>
<p>Once we understand our data, we will demonstrate how <strong>key features of the HaplotypeCaller</strong> enable it to produce better results than position-based callers like UnifiedGenotyper. In particular, we’ll show how <strong>local assembly of haplotypes and realignment of reads</strong> are crucial to producing superior indel calls. Along the way we’ll show you useful tips and tricks for <strong>troubleshooting variant calls</strong> with HaplotypeCaller and the IGV genome browser.</p>
<p>All this will build up to demonstrating the <strong>GVCF workflow for joint variant analysis</strong>, as applied to a trio of whole-genome samples. We hope to convince you that this workflow has substantial practical advantages over a joint analysis that is achieved by calling variants simultaneously on all samples, while producing <strong>results that are just as good</strong> or even better.</p>
<p>The tutorial dataset is available for public download <a href="https://drive.google.com/folderview?id=0BwTg3aXzGxEDNTF3M2hhSnBPU2s&amp;usp=sharing">here</a>.</p>
<hr />
<h3>Table of Contents</h3>
<ol>
<li>WORKING WITH DATASETS FROM DIFFERENT EXPERIMENTAL DESIGNS
1.1 <a href="#1.1">The genome reference: b37</a>
1.2 <a href="#1.2">The test sample: NA12878 Whole-Genome Sequence (WGS)</a>
1.3 <a href="#1.3">For comparison: NA12878 Exome Sequence</a>
1.4 <a href="#1.4">Another comparison: NA12878 RNAseq </a></li>
<li>DIAGNOSING UNKNOWN BAMS
2.1 <a href="#2.1">View header and check read groups</a>
2.2 <a href="#2.2">Validate the file</a></li>
<li>VARIANT DISCOVERY
3.1 <a href="#3.1">Call variants with a position-based caller: UnifiedGenotyper</a>
3.2 <a href="#3.2">Call variants with HaplotypeCaller</a>
&emsp;3.2.1 <a href="#3.2.1">View realigned reads and assembled haplotypes</a>
&emsp;3.2.2 <a href="#3.2.2">Run more samples</a>
3.3 <a href="#3.3">Run HaplotypeCaller on a single bam file in GVCF mode</a>
&emsp;3.3.1 <a href="#3.3.1">View resulting GVCF file in the terminal</a>
&emsp;3.3.2 <a href="#3.3.2">View variants in IGV</a><a name="1.1"></a>
&emsp;3.3.3 <a href="#3.3.3">Run joint genotyping on the CEU Trio GVCFs to generate the final VCF</a>
&emsp;3.3.4 <a href="#3.3.4">View variants in IGV and compare callsets</a></li>
</ol>
<hr />
<h3>1 WORKING WITH DATASETS FROM DIFFERENT EXPERIMENTAL DESIGNS</h3>
<h3>1.1 The genome reference: b37</h3>
<p>We are using a version of the b37 human genome reference containing only a subset of chromosome 20, which we prepared specially for this tutorial in order to provide a reasonable bundle size for download. It is accompanied by its index and sequence dictionary.</p>
<table class="table table-striped">
<thead>
<tr>
<th></th>
<th></th>
</tr>
</thead>
<tbody>
<tr>
<td><font face="Courier New" size="2">ref/</font></td>
<td></td>
<td></td>
</tr>
<tr>
<td><font face="Courier New" size="2">human_g1k_b37_20.fasta</font></td>
<td>&emsp;</td>
<td><font face="Courier New" size="2">genome reference</font></td>
</tr>
<tr>
<td><font face="Courier New" size="2">human_g1k_b37_20.fasta.fai</font></td>
<td>&emsp;</td>
<td><font face="Courier New" size="2">fasta index</font></td>
</tr>
<tr>
<td><font face="Courier New" size="2">human_g1k_b37_20.dict</font></td>
<td>&emsp;</td>
<td><font face="Courier New" size="2">sequence dictionary</font></td>
</tr>
</tbody>
</table>
<p>Open up IGV, and load the <strong>Human (1kg, b37+decoy)</strong> reference available on the IGV server (Genomes&gt;Load Genome from Server). We use this reference in IGV because it has a pre-loaded gene track, whereas our custom chromosome-20-only reference does not.</p>
<p><a name="1.2"></a><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/50/ea77178d2c4407230d0a0282777951.png" /></p>
<h3>1.2 The test sample: NA12878 Whole-Genome Sequence (WGS)</h3>
<p>The biological sample from which the example sequence data was obtained comes from individual NA12878, a member of a 17 sample collection known as CEPH Pedigree 1463, taken from a family in Utah, USA. A trio of two parents and one child from this data set is often referred to as the CEU Trio and is widely used as an evaluation standard (e.g. in the Illumina Platinum Genomes dataset). Note that an alternative trio constituted of the mother (NA12878) and her parents is often also referred to as a CEU Trio. Our trio corresponds to the 2nd generation and one of the 11 grandchildren. </p>
<p>We will begin with a bit of data exploration by looking at the following BAM files derived from NA12878: </p>
<ol>
<li>
<p><code>NA12878_wgs_20.bam</code></p>
<p>Whole genome sequence (WGS) dataset, paired-end 151 bp reads sequenced on Illumina HiSeqX and fully pre-processed according to the GATK Best Practices for germline DNA.</p>
</li>
<li>
<p><code>NA12878_rnaseq_20.bam</code></p>
<p>RNAseq dataset, paired-end 75 bp reads sequenced on Illumina HiSeqX and aligned using STAR 2-pass according to the GATK Best Practices for RNAseq. </p>
</li>
<li>
<p><code>NA12878_ICE_20.bam</code></p>
<p>Exome dataset, Illumina Capture Exome (ICE) library, paired-end 76 bp reads sequenced on Illumina HiSeqX, fully pre-processed according to the GATK Best Practices for germline DNA.</p>
</li>
<li>
<p><code>NA12878_NEX_20.bam</code></p>
<p>Exome dataset, Illumina Nextera Rapid Capture Exome (NEX) library, paired-end 76 bp reads sequenced on Illumina HiSeqX, fully pre-processed according to the GATK Best Practices for germline DNA.</p>
</li>
</ol>
<p>The sequence data files have been specially prepared as well to match our custom chromosome 20-only reference. They only contain data on chromosome 20, in two pre-determined intervals of interest ranging from positions 20:10,000,000-10,200,000 and 20:15,800,000-16,100,00 to keep file sizes down. </p>
<p>Let’s start by loading the DNA WGS sample of NA12878 (<code>bams/exp_design/NA12878_wgs_20.bam</code>), as shown in the screenshots below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/10/48eeee26b88c6acafde8a76e8da4dd.png" />
<p>Initially you will not see any data displayed. You need to zoom in to a smaller region for IGV to start displaying reads. You can do that by using the -/+ zoom controls, or by typing in some genome regions coordinates. Here, we’ll zoom into a predetermined interval of interest, so type <strong>20:16,029,744-16,030,079</strong> into the coordinates box. Once you hit the <code>[Go]</code> button, you should see something like this:</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ee/bfd3a3479f37abd2d6e421ec67aa74.png" />
<p>The top track shows depth of coverage, i.e. the amount of sequence reads present at each position. The mostly grey horizontal bars filling the viewport are the reads. Grey means that those bases match the reference, while colored stripes or base letters (depending on your zoom level) indicate mismatches. You will also see some reads with mapping insertions and deletions, indicated by purple <code>I</code> symbols and crossed-out gaps, respectively.</p>
<p><a name="1.3"></a></p>
<blockquote>
<p><b>TOOL TIP</b>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/19/e73d751f1284ba7c8ea6434f80b12c.png" width="425" align=left hspace="10" /> Read details are shown when you hover over them with your mouse--which can be convenient when troubleshooting, but gets annoying quickly. To turn it off, Click the yellow speech bubble in the toolbar and select “Show details on click”.</p>
</blockquote>
<h3>1.3 For comparison: NA12878 Exome Sequence</h3>
<p>Next, let’s load our two Exome data sets (File&gt;Load from File), <code>NA12878_ICE_20.bam</code> and <code>NA12878_NEX_20.bam</code>, and go to position <strong>20:15,873,697-15,875,416</strong>.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/f2/420a54e32d7e5c9f1cbdda8850b513.png" />
<p>You can see from the coverage graph that the ICE sample has more breadth and depth of coverage at this target site, in comparison to the NEX sample. This directly affects our ability to call variants in the leftmost peak, since ICE provides much more depth and NEX has a particularly lopsided distribution of coverage at that site. That’s not to say that ICE is better in general--just that for this target site, in this <a name="1.4"></a>sequencing run, it provided more even coverage. The overarching point here is that exome kits are not all equivalent and you should evaluate which kit provides the results you need in the regions you care about, before committing to a particular kit for a whole project. As a corollary, comparing exome datasets generated with different kits can be complicated and requires careful evaluation.</p>
<h3>1.4 Another comparison: NA12878 RNAseq</h3>
<p>Lastly, let’s load (File&gt;Load from File) the aligned RNAseq dataset that we have for NA12878 (<code>NA12878_rnaseq_20.bam</code>).</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/cb/abc1c5ef312345bd5e83beb39ba498.png" />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/61/a374bd03d52bdf5ce094a3e211425c.png" align=right />
<p><a name="2.1"></a>You’ll notice pale blue lines to the right of center instead of reads. This is because it’s an intronic region! The blue lines connect to reads that are located in the exon. Click on one to see the N operator in the CIGAR string: in the example here, 32M91225N43M indicates that the read covers a 91225 bp intron.</p>
<hr />
<h3>2 DIAGNOSING UNKNOWN BAMS</h3>
<h3>2.1 View header and check read groups</h3>
<p>Now let’s say that you have been brought on to a new project: you will be analyzing sequenced genomes for particular variants in chromosome 20--since you are the chromosome 20 specialist. Your coworker has given you some files that they sequenced a while back. Unfortunately, their lab notebook mostly illegible and lacking in detail where you can read it. So how do you know what’s been done to these files already? Or even if they are good to use still?</p>
<p>Enter Samtools. You can use this tool to open up the bam file your coworker gave you, and check the bam’s record log. Open up your terminal and execute the following:</p>
<pre><code>samtools view -H bams/exp_design/NA12878_wgs_20.bam | grep ‘@RG’</code></pre>
<p>The bam records log information in the header, so we use <code>view -H</code> to ask it to just show us the header. Since we want to see what this sample is, we will also add <code>| grep ‘@RG’</code>, which will only grab the line of the header that starts with @RG.</p>
<blockquote>
<p><font face="Courier New" size ="2">@RG ID:H0164.2  PL:illumina PU:H0164ALXX140820.2    LB:Solexa-272222    PI:0    DT:2014-08-20T00:00:00-0400 SM:NA12878  CN:BI</font></p>
</blockquote>
<p>You can use the read group information to confirm that this file is what your coworker’s notebook scribbles say it is. You can see that it is indeed the NA12878 sample (SM), and the read group ID H0164.2 (ID) matches, etc. After checking that these identifiers match what you can decipher from your coworker’s writing, call Samtools again. This time we will look at <code>@PG</code> to see what tools have been used on this bam file.</p>
<pre><code>samtools view -H bams/exp_design/NA12878_wgs_20.bam | grep ‘@PG’</code></pre>
<p>Again, this only grabs <code>@PG</code> lines from the header, but you will still get a rather long print out in the terminal; we show a single <code>@PG</code> entry below.</p>
<blockquote>
<p><font face="Courier New" size ="2">@PG ID:bwamem   PN:bwamem   VN:0.7.7-r441   CL:/seq/software/picard/1.750/3rd_party/bwa_mem/bwa mem -M -t 10 -p /ref/b37.fasta /dev/stdin  &gt;  /dev/stdout</font></p>
</blockquote>
<p>At the very beginning of each <code>@PG</code> entry, there will be a program ID. From this entry, you can see that <strong>BWA MEM</strong> was run on the bam file your coworker gave you--the rest of the entry describes the specific parameters that the tool was run with. Scanning through all the entries, you should see that your coworker ran GATK IndelRealigner, GATK PrintReads, MarkDuplicates, and BWA MEM. These tools correlate with the pre-processing steps that your coworker told you they took: mapping with BWA MEM, duplicate marking with MarkDuplicates, indel realignment with IndelRealigner, and lastly, BQSR with PrintReads<em>.
<a name="2.2"></a>
<small></em>How does BQSR correspond to PrintReads? Well, PrintReads is the tool used after running BQSR to apply the recalibration to the bam file itself. Since running BaseRecalibrator didn’t modify the bam file, it isn’t recorded in the bam header, but you can infer that it was run because PrintReads shows up in the header.</small></p>
<h3>2.2 Validate the file</h3>
<p>Now satisfied that the file your coworker gave you is properly pre-processed from looking at its header, you want to make sure that the body of the bam file wasn’t broken at some point. We will try <a href="https://www.broadinstitute.org/gatk/blog?id=7567">diagnosing possible problems</a> in the bam using ValidateSamFile.</p>
<pre><code>java -jar picard.jar ValidateSamFile \
    I=input.bam \
    MODE=SUMMARY</code></pre>
<p>Since we don’t know what kind of errors or warnings we will find, we first run the tool in <code>SUMMARY</code> mode. This will output a histogram listing all the errors and warnings in our file.</p>
<blockquote>
<p><font face="Courier New" size ="2">## HISTOGRAM    java.lang.String
Error Type Count
ERROR:MATE_NOT_FOUND   77</font></p>
</blockquote>
<p>That many errors? The file could be badly damaged, but let’s take a closer look. The error here is a MATE_NOT_FOUND, indicating that a read was marked as paired, but that its mate is not found in the file. Now, usually this would be a point of concern, but your coworker told you that this file was subset to a small part of chromosome 20, so it would make sense that some reads mapped within this region and their mates mapped outside the region. </p>
<p><a name="3.1"></a>We can safely ignore this warning. For more details on errors and warnings that ValidateSamFile can produce (since you won’t just be running your coworker’s samples forever), check out <a href="https://www.broadinstitute.org/gatk/guide/article?id=7571">this article</a>. For your coworker’s file, though, you are finally ready to move on to…</p>
<hr />
<h3>3 VARIANT DISCOVERY</h3>
<h3>3.1 Call variants with a position-based caller: UnifiedGenotyper</h3>
<p>You found a (typed!) copy of your coworker's variant discovery protocol, so you want to run their bam file following it. It tells you to run the following command: </p>
<pre><code>java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_UG_calls.vcf \
    -glm BOTH \
    -L 20:10,000,000-10,200,000</code></pre>
<p>Reading from the protocol, you see that <code>-glm BOTH</code> tells the tool to call both indels and SNPs, while <code>-L</code> gives the interval that the bam was subset to--no use wasting time trying to run on the whole genome when you only have data for a small amount.</p>
<p>When the results return, load the original bam file (<code>bams/exp_design/NA12878_wgs_20.bam</code>) and the output VCF (<code>sandbox/NA12878_wgs_20_UG_calls.vcf</code>) in IGV. Zooming to the coordinates <strong>20:10,002,371-10,002,546</strong>, you will see something like the screenshot below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/d3/6745bf8fbf068735411c14b4de40e9.png" />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/36/fda3f5b642938e99a25a5c504b7711.png" align=left width="300" hspace="10" vspace="10"/>
<p>The variant track shows only variant calls--so at this particular site, there is a homozygous SNP call. (You can click on the variant call for more information on it, too.) The bam track below shows the supporting read data that led to a variant call at that site. </p>
<p>Since this laptop screen is so tiny (our budget went to reagents rather than monitors…) and we can’t zoom out any more vertically, right-click on the bam track and select “Collapsed” view.</p>
<p><a name="3.2"></a>This gives us a better overview of what the data looks like in this region: good even coverage, not too much noise in the region, and reasonable allele balance (mostly variant supports the homozygous variant call). Based on the information we see here, this should be a clear variant site.</p>
<h3>3.2 Call variants with HaplotypeCaller</h3>
<p>While preparing for this project, though, you recall hearing about another variant caller: HaplotypeCaller. And, looking on GATK’s website, you see that it recommends calling your variants using HaplotypeCaller over the old UnifiedGenotyper. The new algorithm calls both SNP and indel variants simultaneously via local de-novo assembly of haplotypes in an active region. Essentially, when this variant caller finds a region with signs of variation, it tosses out the old alignment information (from BWA MEM) and performs a local realignment of reads in that region. This makes HaplotypeCaller more accurate in regions that are traditionally difficult to call--such as areas that contain different types of variants close together. Position-based callers like UnifiedGenotyper simply can’t compete.</p>
<p>You decide to re-run your sample with the new variant caller to see if it makes a difference. Tool documentation on the website gives you a basic command to run, and you add your coworker’s interval trick (<code>-L</code>) in as well.</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_calls.vcf \
    -L 20:10,000,000-10,200,000</code></pre>
<p>Load the output VCF (<code>sandbox/NA12878_wgs_20_HC_calls.vcf</code>) in IGV to compare the HC calls to the previously-loaded UG calls.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/3e/f3ea7dc7def23f0a76ef4866bf3f15.png" />
<p>We see that HC called the same C/T SNP as UG, but it also called another variant, a homozygous variant insertion of three T bases. How is this possible when so few reads seem to support an insertion at this position?</p>
<blockquote>
<p><b>TOOL TIP</b>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/29/654cffb26222c4bbacf32e415fed54.png" width="475" align=left hspace="10" /> When you encounter indel-related weirdness, turn on the display of soft-clips, which IGV turns off by default. Go to View &gt; Preferences &gt; Alignments and select “Show soft-clipped bases” </p>
</blockquote>
<p>With soft clip display turned on, the region lights up with variants. This tells us that the aligner (here, BWA MEM) had a lot of trouble mapping reads in the region. It suggests that HaplotypeCaller may have found a different alignment after performing its local graph assembly step. This reassembled region provided HaplotypeCaller with enough support to call the indel that UnifiedGenotyper missed.
<a name="3.2.1"></a></p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/1b/95a1166568fab9e0edc99a7fa74144.png" />
<h3><em>3.2.1 View realigned reads and assembled haplotypes</em></h3>
<p>But we’re not satisfied with “probably” here. Let’s take a peek under the hood of HaplotypeCaller. You find that HaplotypeCaller has a parameter called <code>-bamout</code>, which allows you to ask for the realigned version of the bam. That realigned version is what HaplotypeCaller uses to make its variant calls, so you will be able to see if a realignment fixed the messy region in the original bam.</p>
<p>You decide to run the following command:</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_calls_debug.vcf \
    -bamout sandbox/NA12878_wgs_20.HC_out.bam \
    -forceActive -disableOptimizations \
    -L 20:10,002,371-10,002,546 -ip 100</code></pre>
<p>Since you are only interested in looking at that messy region, you decide to give the tool a narrowed interval with <code>-L 20:10,002,371-10,002,546</code>, with a 100 bp padding on either side using <code>-ip 100</code>. To make sure the tool does perform the reassembly in that region, you add in the <code>-forceActive</code> and <code>-disableOptimizations</code> arguments.</p>
<p>Load the output BAM (<code>sandbox/NA12878_wgs_20.HC_out.bam</code>) in IGV, and switch to Collapsed view once again. You should still be zoomed in on coordinates <strong>20:10,002,371-10,002,546</strong>, and have the original bam track loaded for comparison.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/e3/ba5887413f0f13b527ffd92a1a89b0.png" />
<p>After realignment by HaplotypeCaller (the bottom track), almost all the reads show the insertion, and the messy soft clips from the original bam are gone. Expand the reads in the output BAM (right click&gt;Expanded view), and you can see that all the insertions are in phase with the C/T SNP. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/2b/0f600e7939827adbe896464a7b5a06.png" />
<p>There is more to a BAM than meets the eye--or at least, what you can see in this view of IGV. Right-click on the reads to bring up the view options menu. Select <strong>Color alignments by</strong>, and choose <strong>read group</strong>. Your gray reads should now be colored similar to the screenshot below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/46/4c2802a5aff4d18baac7af7f144771.png" />
<p>Some of the first reads, shown in red at the top of the pile, are not real reads. These represent artificial haplotypes that were constructed by HaplotypeCaller, and are tagged with a special read group identifier, “ArtificialHaplotype,” so they can be visualized in IGV. You can click on an artificial read to see this tag under <strong>RG</strong>.</p>
<p>We see that HaplotypeCaller considered six possible haplotypes, because there is more than one variant in the same ActiveRegion. Zoom out further , and we can see that two ActiveRegions were examined within the scope of the interval we provided (with padding).
<a name="3.2.2"></a></p>
<center><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ee/19be1d025100c211f410161ea1917f.png" width="550" /></center>
<h3><em>3.2.2 Run more samples</em></h3>
<p>You’ve decided that perhaps HaplotypeCaller will work better for your project. However, since you have been working on this protocol update, your coworker found two more samples--they were in a different folder on their computer for reasons you can’t figure out. Regardless, you now need to joint call all the samples together. So, using the same command as before, you’ve tacked on the two additional bam files.</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -I bams/trio-calling/NA12877_wgs_20.bam \
    -I bams/trio-calling/NA12882_wgs_20.bam \
    -o sandbox/NA12878_wgs_20_HC_jointcalls.vcf \
    -L 20:10,000,000-10,200,000</code></pre>
<p><a name="3.3"></a>
You notice after entering that, that HaplotypeCaller takes a much longer time to return than other tasks we have run so far. You decide to check the results of this command later, and do some digging on how to make things go faster.</p>
<h3>3.3 Run HaplotypeCaller on a single bam file in GVCF mode</h3>
<p>Every time your coworker finds a new folder of samples, you’ll have to re-run all the samples using this increasingly slower HaplotypeCaller command. You’ve also been approved for a grant and intend to send your own samples out for sequencing, so there are those to add in as well. You could just wait until you have all the samples gathered, but that could be a while and your PI wants to see some preliminary results soon. You read about a new GATK workflow that lets you make everyone happy: the GVCF workflow. </p>
<p>The first step in variant discovery is to run HaplotypeCaller in GVCF mode on each individual bam file. This is basically running HaplotypeCaller as you did before, but with <code>-ERC GVCF</code> added to the command. You first want to  run HaplotypeCaller in GVCF mode on the NA12878 bam. (In the interest of time, we have supplied the other sample GVCFs in the bundle, but normally you would run them individually in the same way as the first.) This will produce a GVCF file that contains genotype likelihoods for each variant position as well as blocks for each interval where no variant is likely. You’ll see what this looks like more in a minute.
<a name="3.3.1"></a></p>
<pre><code>java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R ref/human_g1k_b37_20.fasta \
    -I bams/exp_design/NA12878_wgs_20.bam \
    -o sandbox/NA12878_wgs_20.g.vcf \
    -ERC GVCF \
    -L 20:10,000,000-10,200,000</code></pre>
<h3><em>3.3.1 View resulting GVCF file in the terminal</em></h3>
<p>Since a GVCF is a new file type for your workflow, let’s take a look at the actual content first. You can do this in the terminal by typing this command:</p>
<p>more sandbox/NA12878_wgs_20.g.vcf</p>
<p>As you scroll through the file (hit <code>[ENTER]</code> to scroll, <code>[CTRL]</code>+<code>[C]</code> to exit), note the NON_REF allele defined in the header.</p>
<blockquote>
<p><font face="Courier New" size ="2">##ALT=&lt;ID=NON_REF,Description=”Represents any possible alternative allele at this location”&gt;</font></p>
</blockquote>
<p>Also note the GVCF blocks defined later in the header. The reference (non-variant) blocks are recorded in the GVCF file, in blocks separated by genotype quality.</p>
<blockquote>
<p><font face="Courier New" size ="2">##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)</font>
<font face="Courier New" size ="2">##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)</font>
<font face="Courier New" size ="2">##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)</font>
<font face="Courier New" size ="2">##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)</font></p>
</blockquote>
<p>Finally, while scrolling through the records, we can see the <font face="Courier New" size ="2">reference blocks</font> and <font face="Courier New" size ="2"><b>variant sites</b></font>.</p>
<blockquote>
<p><font face="Courier New" size ="2">20  10000115    .   G   <NON_REF>   .   .   END=10000116    GT:DP:GQ:MIN_DP:PL  0/0:25:69:25:0,69,1035
<b>20  10000117    .   C   T,<NON_REF> 262.77  .   BaseQRankSum=-0.831;ClippingRankSum=-0.092;DP=23;MLEAC=1,0;MLEAF=0.500,0.00;MQ=60.47;MQRankSum=1.446;ReadPosRankSum=0.462   GT:AD:DP:GQ:PL:SB   0/1:11,12,0:23:99:291,0,292,324,327,652:9,2,9,3</b>
20 10000118    .   T   <NON_REF>   .   .   END=10000123    GT:DP:GQ:MIN_DP:PL  0/0:25:63:24:0,63,945</font></p>
</blockquote>
<p><a name="3.3.2"></a>Every site in the interval we analyzed is represented here--whether it be by a variant call, a reference call, or a reference block. This helps to distinguish between a “no call” (we don’t have enough data to make a call) and a “reference call” (we have evidence that the site matches the reference).</p>
<h3><em>3.3.2 View variants in IGV</em></h3>
<p>Now, text in a terminal window can be rather hard to read, so let’s take a look at the GVCFs in IGV. Start a new session to clear your IGV screen, then load the three GVCFs (<code>sandbox/NA12878_wgs_20.g.vcf</code>, <code>gvcfs/NA12877_wgs_20.g.vcf</code>, <code>gvcfs/NA12882_wgs_20.g.vcf</code>). You should already be zoomed in on <strong>20:10,002,371-10,002,546</strong> from our previous section, and see this:</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/49/98eaa9c2db94b98f9f148948bd21b4.png" />
<p><a name="3.3.3"></a>Notice anything different from the VCF? Along with the colorful variant sites, you see many gray blocks in the GVCF representing the non-variant intervals. Most of the gray blocks are next to each other, but are not grouped together, because they belong to different GQ blocks. The chief difference between the GVCF here and the next step’s VCF is the lack of reference blocks (the gray bits). Only very low-confidence variant sites will be removed in the VCF, based on the QUAL score.</p>
<h3><em>3.3.3 Run joint genotyping on the CEU Trio GVCFs to generate the final VCF</em></h3>
<p>The last step is to joint call all your GVCF files using the GATK tool GenotypeGVCFs. After looking in the tool documentation, you run this command:</p>
<pre><code>java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \
    -R ref/human_g1k_b37_20.fasta \
    -V sandbox/NA12878_wgs_20.g.vcf \
    -V gvcfs/NA12877_wgs_20.g.vcf \
    -V gvcfs/NA12882_wgs_20.g.vcf \
    -o sandbox/CEUTrio_wgs_20_GGVCFs_jointcalls.vcf \
    -L 20:10,000,000-10,200,000</code></pre>
<p><a name="3.3.4"></a>
That returned much faster than the HaplotypeCaller step--and a good thing, too, since this step is the one you’ll need to re-run every time your coworker finds a “new” sample buried in their messy file structure. But does calling this way really give you good results? Let’s take a look.</p>
<h3><em>3.3.4 View variants in IGV and compare callsets</em></h3>
<p>Load the joint called VCF from normal HaplotypeCaller, section 3.2.1 (<code>sandbox/NA12878_wgs_20_HC_jointcalls.vcf</code>), and GenotypeGVCFs, section 3.3.3 (<code>sandbox/CEUTrio_wgs_20_GGVCFs_jointcalls.vcf</code>). Change your view to look at <strong>20:10,002,584-10,002,665</strong>, and you will see:</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/be/01d8f53fcd5e27bcfef44366cafed8.png" />
<p>At this site, the father NA12877 is heterozygous for a G/T SNP, and the mother, NA12878, and son, NA12882, are homozygous variant for the same SNP. These calls match up, and you figure that the calls between GenotypeGVCFs and HaplotypeCaller, when run in multisample mode, are essentially equivalent. (And if you did some digging, you would find some marginal differences in borderline calls.) However, the GVCF workflow allows you to be more flexible. Every time your PI wants an update on the project, you can simply re-run the quick GenotypeGVCFs step on all the samples you have gathered so far. The expensive and time-consuming part of calculating genotype likelihoods only needs to be done once on each sample, so you won’t have to spend all your grant money on compute to rerun the whole cohort every time you have a new sample.</p>
<p>You have successfully run your coworker’s samples, and you’ve found that the most effective workflow for you is the most recent GVCF workflow. Your next step takes you to filtering the callset with either VQSR or hard filters--but you decide to take a break before tackling the next part of the workflow.</p>