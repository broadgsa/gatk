## (howto) Recalibrate base quality scores = run BQSR

http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

<h4>Objective</h4>
<p>Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.</p>
<h4>Prerequisites</h4>
<ul>
<li>TBD</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Analyze patterns of covariation in the sequence dataset </li>
<li>Do a second pass to analyze covariation remaining after recalibration </li>
<li>Generate before/after plots</li>
<li>Apply the recalibration to your sequence data</li>
</ol>
<hr />
<h3>1. Analyze patterns of covariation in the sequence dataset</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I input_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -o recal_data.table </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a GATKReport file called <code>recal_data.table</code> containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data. </p>
<p>It is imperative that you provide the program with a set of known sites, otherwise it will refuse to run. The known sites are used to build the covariation model and estimate empirical base qualities. For details on what to do if there are no known sites available for your organism of study, please see the online GATK documentation. </p>
<p>Note that <code>-L 20</code> is used here and in the next steps to restrict analysis to only chromosome 20 in the b37 human genome reference build. To run against a different reference, you may need to change the name of the contig according to the nomenclature used in your reference. </p>
<hr />
<h3>2. Do a second pass to analyze covariation remaining after recalibration</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I realigned_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -BQSR recal_data.table \ 
    -o post_recal_data.table </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates another GATKReport file, which we will use in the next step to generate plots. Note the use of the <code>-BQSR</code> flag, which tells the GATK engine to perform on-the-fly recalibration based on the first recalibration data table. </p>
<hr />
<h3>3. Generate before/after plots</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T AnalyzeCovariates \ 
    -R reference.fa \ 
    -L 20 \ 
    -before recal_data.table \
    -after post_recal_data.table \
    -plots recalibration_plots.pdf</code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This generates a document called <code>recalibration_plots.pdf</code> containing plots that show how the reported base qualities match up to the empirical qualities calculated by the BaseRecalibrator. Comparing the <strong>before</strong> and <strong>after</strong> plots allows you to check the effect of the base recalibration process before you actually apply the recalibration to your sequence data. For details on how to interpret the base recalibration plots, please see the online GATK documentation. </p>
<hr />
<h3>4. Apply the recalibration to your sequence data</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T PrintReads \ 
    -R reference.fa \ 
    -I input_reads.bam \ 
    -L 20 \ 
    -BQSR recal_data.table \ 
    -o recal_reads.bam </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a file called <code>recal_reads.bam</code> containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag <code>–emit_original_quals</code> to the PrintReads command, in which case the original qualities will also be written in the file, tagged <code>OQ</code>.</p>
<p>Notice how this step uses a very simple tool, PrintReads, to apply the recalibration. What’s happening here is that we are loading in the original sequence data, having the GATK engine recalibrate the base qualities on-the-fly thanks to the <code>-BQSR</code> flag (as explained earlier), and just using PrintReads to write out the resulting data to the new file.</p>