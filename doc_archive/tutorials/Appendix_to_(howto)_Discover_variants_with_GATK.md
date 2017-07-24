## Appendix to (howto) Discover variants with GATK

http://gatkforums.broadinstitute.org/gatk/discussion/7870/appendix-to-howto-discover-variants-with-gatk

<h2>GATK TUTORIAL :: Variant Discovery :: Appendix</h2>
<p><strong>June 2016 - GATK 3.6</strong></p>
<p>This document is an appendix to the <a href="https://www.broadinstitute.org/gatk/guide/article?id=7869">GATK Tutorial :: Variant Discovery module worksheet</a>. It contains a summary introduction to the scientific context of the tutorial. </p>
<hr />
<h3>Table of Contents</h3>
<ol>
<li><a href="#1">GATK BEST PRACTICES</a></li>
<li><a href="#2">WHAT IS JOINT ANALYSIS?</a></li>
<li><a href="#3">FLAWS OF JOINT ANALYSIS</a>
3.1 <a href="#3.1">The N+1 problem</a><a name="1"></a>
3.2 <a href="#3.2">Really bad scaling</a></li>
<li><a href="#4">THE GVCF WORKFLOW</a></li>
</ol>
<hr />
<h3>1 GATK BEST PRACTICES</h3>
<p>The GATK Best Practices workflows provide step-by-step recommendations for performing variant discovery analysis in high-throughput sequencing (HTS) data. The following diagram illustrates the GATK Best Practices workflow for germline SNP and Indel discovery in whole genomes and exomes. It includes three phases: pre-processing, variant discovery, and callset refinement.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/00/72ce1625d9d0cc4d35445a0aaf34ea.png" />
<p><strong>Figure 1: Best Practices workflow for germline SNP and Indel discovery in whole-genomes and exomes.</strong></p>
<p><strong>Pre-Processing</strong> starts from raw sequence data, either in FASTQ or uBAM format, and produces analysis-ready BAM files. Processing steps include alignment to a reference genome as well as some data cleanup operations to correct for technical biases and make the data suitable for analysis.</p>
<p><strong>Variant Discovery</strong> starts from analysis-ready BAM files and produces a callset in VCF format. Processing involves identifying sites where one or more individuals display possible genomic variation, and applying filtering methods appropriate to the experimental design. The <strong>Best Practices version 3.x</strong> include key innovations that enable <strong>joint analysis of multiple samples</strong> in a way that is <strong>scalable</strong> and allows <strong>incremental processing</strong> of the sequencing data. Those innovations are the focus of this tutorial.</p>
<p><strong>Callset Refinement</strong> starts and ends with a VCF callset. Processing involves using metadata such as previously validated callsets to assess and improve genotyping accuracy, attach additional information and evaluate the overall quality of the callset.
<a name="2"></a>
Learn more about the GATK Best Practices <a href="https://www.broadinstitute.org/gatk/guide/best-practices">here</a>.</p>
<hr />
<h3>2 WHAT IS JOINT ANALYSIS?</h3>
<p>In this context, joint analysis means that we consider evidence from multiple samples in order to determine the genotype of each sample at each site, rather than looking at only one sample at a time in isolation. Considering evidence from multiple samples empowers variant discovery and allows us to detect variants with great sensitivity and genotype samples as accurately as possible. Specifically, we have determined that joint analysis conveys the following benefits:</p>
<ul>
<li>Clearer distinction between homozygous reference sites and sites with missing data</li>
<li>Greater sensitivity for low-frequency variants in low-coverage data</li>
<li>Greater ability to filter out false positives without losing sensitivity</li>
</ul>
<p>There are specific data contexts in which performing joint analysis makes an especially important difference. Two such cases are illustrated below.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/f1/2a748628e91c91ec9f2cbdfab70a74.png" />
<p><strong>Figure 2: Two cases where joint analysis provides important information that improves either the genotype determination or the interpretation of results.</strong></p>
<p><strong>Left: Power of joint analysis in finding mutations at low coverage sites.</strong> The variant allele is present in only two of the N samples, in both cases with such low coverage that the variant is not callable when processed separately. Joint calling allows evidence to be accumulated over all samples and renders the variant callable. </p>
<p><strong>Right: Importance of joint analysis to square off the genotype matrix, using an example of two disease-relevant variants.</strong> If we call these samples independently and produce a variants-only output, neither sample will have records for these two sites, for different reasons: the first sample is homozygous reference while the second sample has no data. Therefore, <a name="3"></a>merging the results from single sample calling will incorrectly treat both of these samples identically as being non-informative.</p>
<p>Learn more about joint analysis <a href="https://www.broadinstitute.org/gatk/guide/article?id=4150">here</a>.</p>
<hr />
<h3>3 FLAWS OF JOINT ANALYSIS</h3>
<img src= "https://us.v-cdn.net/5019796/uploads/FileUpload/a9/1b2ed7501e3a82ee3f6770dbfb89a4.png" align=right width=225 />
<p><a name="3.1"></a>Traditionally, joint analysis was achieved by calling variants jointly across all sample BAMs at the same time, generating a single call set for the entire cohort in a single step. </p>
<p>However, that method suffers from two major flaws: the <strong>N+1 problem</strong> and <strong>really bad scaling</strong>.</p>
<h3>3.1 The N+1 problem</h3>
<p>When you’re getting a large-ish number of samples sequenced (especially clinical samples), you typically get them in small batches over an extended period of time. In the past, this was handled by doing batch calling, i.e. analyze the samples in batches and combine the resulting VCF callsets as they become available. But that’s not a true joint analysis, and it doesn’t give you the same significant gains that calling variants jointly can yield (on top of producing batch effects). If <a name="3.2"></a>you wanted to do a true joint analysis using the multisample variant calling approach, you have to re-call all samples from scratch every time you get even one new sample sequence. And the more you add samples, the more computationally intensive it gets, bringing us to the next problem: really bad scaling.</p>
<h3>3.2 Really bad scaling</h3>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/b7/e1f4634cf066cbd389cc0355053ab2.png" align=right width=450 /> Calling variants jointly across samples scales very badly. This is because the calculations involved in variant calling (especially by sophisticated tools like the HaplotypeCaller that perform a graph assembly step) become exponentially more computationally costly as you add samples to the cohort. If you don't have a lot of compute <a name="4"></a>available, you run into limitations very quickly. Even at Broad, where we have fairly ridiculous amounts of compute available, we can't brute-force our way through the numbers for the large cohort sizes that we're called on to handle like the 92,000 exomes of the ExAC dataset (see <a href="http://exac.broadinstitute.org/">this page</a>).</p>
<hr />
<h3>4 THE GVCF WORKFLOW</h3>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/83/fa236e7e23016d0a33fa3e07f1b086.png" align=right width=225 /> The good news is that you don’t actually have to call variants on all your samples together to perform a joint analysis. We have developed a workflow that allows us to <strong>decouple</strong> the initial identification of potential variant sites, i.e. the <strong>variant calling</strong>, from the <strong>genotyping</strong> step, which is the only part that really needs to be done jointly. Since GATK 3.0, you can use the HaplotypeCaller to call variants individually per-sample in a special mode invoked by adding -ERC GVCF to your command line, generating an intermediate file called a <strong>GVCF</strong> (for Genomic VCF). You then run a <strong>joint genotyping</strong> step on all the GVCF files generated for the samples in the cohort. This achieves what we call incremental joint discovery, providing you with all the benefits of classic joint calling (as described below) without the drawbacks.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ec/e76d8cd1decb03d8f8e89c5dba7f7f.png" width=400 align=left hspace="10" />
<p>&#32;
&#32;
&#32;
<strong>Figure 4. The new approach to joint analysis allows incremental processing of samples and scales much better than the traditional approach of calling variants on all samples simultaneously.</strong>
&#32;
&#32;
&#32;
&#32;
&#32;</p>
<hr />
<h3>Conclusion</h3>
<p><strong>This uniquely innovative workflow solves both the scaling problems and the N+1 problem that plague traditional methods of joint analysis.</strong></p>
<p>From here on out we will refer to this <strong>single-sample calling + joint genotyping workflow</strong> as <strong>the GVCF workflow</strong> because it involves the intermediate GVCF file, which uniquely distinguishes it from other methods.</p>