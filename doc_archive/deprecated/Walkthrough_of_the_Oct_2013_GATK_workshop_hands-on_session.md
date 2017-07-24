## Walkthrough of the Oct 2013 GATK workshop hands-on session

http://gatkforums.broadinstitute.org/gatk/discussion/3366/walkthrough-of-the-oct-2013-gatk-workshop-hands-on-session

<h4>Note: the exact data files we used in this tutorial are no longer available. However, you can use the files in the resource bundle to work through this tutorial. You may need to adapt the filenames accordingly.</h4>
<hr />
<h3>Map and mark duplicates</h3>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2799/howto-map-and-mark-duplicates">http://gatkforums.broadinstitute.org/discussion/2799/howto-map-and-mark-duplicates</a></p>
<p><em>Starting with aligned (mapped)  and deduplicated (dedupped) reads in .sam file to save time.</em></p>
<h4>- Generate index</h4>
<p>Create an index file to enable fast seeking through the file.</p>
<pre><code class="pre_md">java -jar BuildBamIndex.jar I= dedupped_20.bam</code class="pre_md"></pre>
<h4>- Prepare reference to work with GATK</h4>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk">http://gatkforums.broadinstitute.org/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk</a></p>
<p>Create a dictionary file and index for the reference.</p>
<pre><code class="pre_md">java -jar CreateSequenceDictionary.jar R=human_b37_20.fasta O=human_b37_20.dict

samtools faidx human_b37_20.fasta </code class="pre_md"></pre>
<hr />
<h3>Getting to know GATK</h3>
<h4>- Run a simple walker: CountReads</h4>
<p>Identify basic syntax, console output: version, command recap line, progress estimates, result if applicable.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T CountReads -R human_b37_20.fasta -I dedupped_20.bam -L 20</code class="pre_md"></pre>
<h4>- Add a filter to count how many duplicates were marked</h4>
<p>Look at filtering summary.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T CountReads -R human_b37_20.fasta -I dedupped_20.bam -L 20 -rf DuplicateRead</code class="pre_md"></pre>
<h4>- Demonstrate how to select a subset of read data</h4>
<p>This can come in handy for bug reports.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T PrintReads -R human_b37_20.fasta -I dedupped_20.bam -L 20:10000000-11000000 -o snippet.bam</code class="pre_md"></pre>
<h4>- Demonstrate the equivalent for variant calls</h4>
<p>Refer to docs for many other capabilities including selecting by sample name, up to complex queries.</p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T SelectVariants -R human_b37_20.fasta -V dbsnp_b37_20.vcf -o snippet.vcf -L 20:10000000-11000000</code class="pre_md"></pre>
<hr />
<h3>Back to data processing</h3>
<h4>- Realign around Indels</h4>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2800/howto-perform-local-realignment-around-indels">http://gatkforums.broadinstitute.org/discussion/2800/howto-perform-local-realignment-around-indels</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R human_b37_20.fasta -I dedupped_20.bam -known indels_b37_20.vcf -o target_intervals.list -L 20 

java -jar GenomeAnalysisTK.jar -T IndelRealigner -R human_b37_20.fasta -I dedupped_20.bam -known indels_b37_20.vcf -targetIntervals target_intervals.list -o realigned_20.bam -L 20 </code class="pre_md"></pre>
<h4>- Base recalibration</h4>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr">http://gatkforums.broadinstitute.org/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human_b37_20.fasta -I realigned_20.bam -knownSites dbsnp_b37_20.vcf -knownSites indels_b37_20.vcf -o recal_20.table -L 20

java -jar GenomeAnalysisTK.jar -T PrintReads -R human_b37_20.fasta -I realigned_20.bam -BQSR recal_20.table -o recal_20.bam -L 20

java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human_b37_20.fasta -I recalibrated_20.bam -knownSites dbsnp_b37_20.vcf -knownSites indels_b37_20.vcf -o post_recal_20.table -L 20

java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R human_b37_20.fasta -before recal_20.table -after post_recal_20.table -plots recalibration_plots.pdf -L 20 </code class="pre_md"></pre>
<h4>- ReduceReads</h4>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2802/howto-compress-read-data-with-reducereads">http://gatkforums.broadinstitute.org/discussion/2802/howto-compress-read-data-with-reducereads</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T ReduceReads -R human_b37_20.fasta -I recalibrated_20.bam -o reduced_20.bam -L 20 </code class="pre_md"></pre>
<h4>- HaplotypeCaller</h4>
<p><a href="http://gatkforums.broadinstitute.org/discussion/2803/howto-call-variants-on-a-diploid-genome-with-the-haplotypecaller">http://gatkforums.broadinstitute.org/discussion/2803/howto-call-variants-on-a-diploid-genome-with-the-haplotypecaller</a></p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_b37_20.fasta -I reduced_20.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o variants_20.vcf -L 20 </code class="pre_md"></pre>