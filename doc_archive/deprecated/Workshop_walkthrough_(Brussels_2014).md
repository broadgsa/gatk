## Workshop walkthrough (Brussels 2014)

http://gatkforums.broadinstitute.org/gatk/discussion/4327/workshop-walkthrough-brussels-2014

<h4>Note: this is a walkthrough of a hands-on GATK tutorial given at the Royal Institute of Natural Sciences on June 26, 2014 in Brussels, Belgium. It is intended to be performed with version 3.1-2 of the GATK and the corresponding data bundle.</h4>
<h3>Data files</h3>
<p>We start with a BAM file called &quot;NA12878.wgs.1lib.bam&quot; (along with its index, &quot;NA12878.wgs.1lib.bai&quot;) containing Illumina sequence reads from our favorite test subject, NA12878, that have been mapped using BWA-mem and processed using Picard tools according to the instructions here:</p>
<p><a href="http://www.broadinstitute.org/gatk/guide/article?id=2799">http://www.broadinstitute.org/gatk/guide/article?id=2799</a></p>
<p>Note that this file only contains sequence for a small region of chromosome 20, in order to minimize the file size and speed up the processing steps, for demonstration purposes. Normally you would run the steps in this tutorial on the entire genome (or exome). </p>
<p>This subsetted file was prepared by extracting read group 20GAV.1 from the CEUTrio.HiSeq.WGS.b37.NA12878.bam that is available in our resource bundle, using the following command:</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T PrintReads -R human_g1k_v37.fasta -I CEUTrio.HiSeq.WGS.b37.NA12878.bam -o NA12878.wgs.1lib.bam -L 20 -rf SingleReadGroup -goodRG 20GAV.1</code class="pre_md"></pre>
<p>(We'll explain later in the tutorial how to use this kind of utility function to manipulate BAM files.)</p>
<p>We also have our human genome reference, called &quot;human_g1k_v37.fasta&quot;, which has been prepared according to the instructions here:</p>
<p><a href="http://www.broadinstitute.org/gatk/guide/article?id=2798">http://www.broadinstitute.org/gatk/guide/article?id=2798</a></p>
<p>We will walk through both of these tutorials to explain the processing, but without actually running the steps to save time.</p>
<p>And finally we have a few resource files containing known variants (dbsnp, mills indels). These files are all available in the resource bundle on our FTP server. See here for access instructions:</p>
<p><a href="http://www.broadinstitute.org/gatk/guide/article?id=1215">http://www.broadinstitute.org/gatk/guide/article?id=1215</a></p>
<hr />
<h2>DAY 1</h2>
<h3>Prelude: BAM manipulation with Picard and Samtools</h3>
<h4>- Viewing a BAM file information</h4>
<p>See also the Samtools docs:</p>
<p><a href="http://samtools.sourceforge.net/samtools.shtml">http://samtools.sourceforge.net/samtools.shtml</a>  </p>
<h4>- Reverting a BAM file</h4>
<p>Clean the BAM we are using of previous GATK processing using this Picard command:</p>
<pre><code class="pre_md">java -jar RevertSam.jar I=NA12878.wgs.1lib.bam O=aligned_reads_20.bam RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=false SORT_ORDER=coordinate</code class="pre_md"></pre>
<p>Note that it is possible to revert the file to FastQ format by setting REMOVE_ALIGNMENT_INFORMATION=true, but this method leads to biases in the alignment process, so if you want to do that, the better method is to follow the instructions given here:</p>
<p><a href="http://www.broadinstitute.org/gatk/guide/article?id=2908">http://www.broadinstitute.org/gatk/guide/article?id=2908</a></p>
<p>See also the Picard docs:</p>
<p><a href="http://picard.sourceforge.net/command-line-overview.shtml">http://picard.sourceforge.net/command-line-overview.shtml</a>  </p>
<h3>Mark Duplicates</h3>
<p>See penultimate step of <a href="http://www.broadinstitute.org/gatk/guide/article?id=2799">http://www.broadinstitute.org/gatk/guide/article?id=2799</a></p>
<p>After a few minutes, the file (which we'll call &quot;dedupped_20.bam&quot;) is ready for use with GATK.</p>
<h3>Interlude: tour of the documentation, website, forum etc. Also show how to access the bundle on the FTP server with FileZilla.</h3>
<h3>Getting to know GATK</h3>
<p>Before starting to run the GATK Best Practices, we are going to learn about the basic syntax of GATK, how the results are output, how to interpret error messages, and so on.</p>
<h4>- Run a simple walker: CountReads</h4>
<p>Identify basic syntax, console output: version, command recap line, progress estimates, result if applicable.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T CountReads -R human_g1k_v37.fasta -I dedupped_20.bam -L 20</code class="pre_md"></pre>
<h4>- Add a filter to count how many duplicates were marked</h4>
<p>Look at the filtering summary.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T CountReads -R human_g1k_v37.fasta -I dedupped_20.bam -L 20 -rf DuplicateRead</code class="pre_md"></pre>
<h4>- Demonstrate how to select a subset of read data</h4>
<p>This can come in handy for bug reports.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T PrintReads -R human_g1k_v37.fasta -I dedupped_20.bam -L 20:10000000-11000000 -o snippet.bam</code class="pre_md"></pre>
<p>Also show how a bug report should be formatted and submitted. See
<a href="http://www.broadinstitute.org/gatk/guide/article?id=1894">http://www.broadinstitute.org/gatk/guide/article?id=1894</a></p>
<h4>- Demonstrate the equivalent for variant calls</h4>
<p>Refer to docs for many other capabilities including selecting by sample name, up to complex queries.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R human_g1k_v37.fasta -V dbsnp_b37_20.vcf -o snippet.vcf -L 20:10000000-11000000</code class="pre_md"></pre>
<p>See <a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_SelectVariants.html">http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_SelectVariants.html</a></p>
<hr />
<h3>GATK Best Practices for data processing (DNA seq)</h3>
<p>These steps should typically be performed per lane of data. Here we are running the tools on a small slice of the data, to save time and disk space, but normally you would run on the entire genome or exome. This is especially important for BQSR, which does not work well on small amounts of data.</p>
<p>Now let's pick up where we left off after Marking Duplicates.</p>
<h4>- Realign around Indels</h4>
<p>See <a href="http://gatkforums.broadinstitute.org/discussion/2800/howto-perform-local-realignment-around-indels">http://gatkforums.broadinstitute.org/discussion/2800/howto-perform-local-realignment-around-indels</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R human_g1k_v37.fasta -I dedupped_20.bam -known Mills_and_1000G_gold_standard.indels.b37 -o target_intervals.list -L 20:10000000-11000000 

java -jar GenomeAnalysisTK.jar -T IndelRealigner -R human_g1k_v37.fasta -I dedupped_20.bam -known Mills_and_1000G_gold_standard.indels.b37.vcf -targetIntervals target_intervals.list -o realigned.bam -L 20:10000000-11000000 </code class="pre_md"></pre>
<h4>- Base recalibration</h4>
<p>See <a href="http://gatkforums.broadinstitute.org/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr">http://gatkforums.broadinstitute.org/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human_g1k_v37.fasta -I realigned_20.bam -knownSites dbsnp_b37_20.vcf -knownSites Mills_and_1000G_gold_standard.indels.b37.vcf -o recal_20.table -L 20:10000000-11000000

java -jar GenomeAnalysisTK.jar -T PrintReads -R human_g1k_v37.fasta -I realigned_20.bam -BQSR recal_20.table -o recal_20.bam -L 20:10000000-11000000

java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human_g1k_v37.fasta -I recalibrated_20.bam -knownSites dbsnp_b37_20.vcf -knownSites Mills_and_1000G_gold_standard.indels.b37.vcf -o post_recal_20.table -L 20:10000000-11000000

java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R human_g1k_v37.fasta -before recal_20.table -after post_recal_20.table -plots recalibration_plots.pdf -L 20:10000000-11000000</code class="pre_md"></pre>
<hr />
<h3>GATK Best Practices for variant calling (DNA seq)</h3>
<h4>- Run HaplotypeCaller in regular mode</h4>
<p>See <a href="http://www.broadinstitute.org/gatk/guide/article?id=2803">http://www.broadinstitute.org/gatk/guide/article?id=2803</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_g1k_v37.fasta -I recal_20.bam -o raw_hc_20.vcf -L 20:10000000-11000000</code class="pre_md"></pre>
<p>Look at VCF in text and in IGV, compare with bam file.</p>
<h4>- Run HaplotypeCaller in GVCF mode (banded and BP_RESOLUTION)</h4>
<p>See <a href="http://www.broadinstitute.org/gatk/guide/article?id=3893">http://www.broadinstitute.org/gatk/guide/article?id=3893</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_g1k_v37.fasta -I recal_20.bam -o raw_hc_20.g.vcf -L 20:10000000-11000000 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000</code class="pre_md"></pre>
<p>Compare to regular VCF.</p>