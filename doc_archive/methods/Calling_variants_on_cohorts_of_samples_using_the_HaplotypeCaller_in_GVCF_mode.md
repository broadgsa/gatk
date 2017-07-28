## Calling variants on cohorts of samples using the HaplotypeCaller in GVCF mode

http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode

<p>This document describes the new approach to joint variant discovery that is available in GATK versions 3.0 and above. For a more detailed discussion of why it's better to perform joint discovery, see this <a href="https://www.broadinstitute.org/gatk/guide/article?id=4150">FAQ article</a>. For more details on how this fits into the overall reads-to-variants analysis workflow, see the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices workflows</a> documentation.</p>
<h3>Overview</h3>
<p>This is the workflow recommended in our Best Practices for performing variant discovery analysis on cohorts of samples.</p>
<p><a href='https://us.v-cdn.net/5019796/uploads/FileUpload/eb/44f317f8850ba74b64ba47b02d1bae.png'><img src='https://us.v-cdn.net/5019796/uploads/FileUpload/eb/44f317f8850ba74b64ba47b02d1bae.png' /></a></p>
<p>In a nutshell, we now call variants individually on each sample using the HaplotypeCaller in <code>-ERC GVCF</code> mode, leveraging the previously introduced <a href="http://www.broadinstitute.org/gatk/guide/article?id=4042">reference model</a> to produce a comprehensive record of genotype likelihoods and annotations for each site in the genome (or exome), in the form of a <a href="http://www.broadinstitute.org/gatk/guide/article?id=4017">gVCF file (genomic VCF)</a>. </p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/9f/f0619642db06b73b599253f42ef2bf.png" />
<p>In a second step, we then perform a joint genotyping analysis of the gVCFs produced for all samples in a cohort.
This allows us to achieve the same results as joint calling in terms of accurate genotyping results, without the computational nightmare of exponential runtimes, and with the added flexibility of being able to re-run the population-level genotyping analysis at any time as the available cohort grows.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a5/ce09f7ae5c11956e2db2f9e763648c.png" />
<p>This is meant to replace the joint discovery workflow that we previously recommended, which involved calling variants jointly on multiple samples, with a much smarter approach that reduces computational burden and solves the &quot;N+1 problem&quot;.</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/6f/aa7af2488c2de8510e556423ee6cfa.png" />
<hr />
<h3>Workflow details</h3>
<p>This is a quick overview of how to apply the workflow in practice. For more details, see the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">Best Practices workflows</a> documentation.</p>
<h4>1. Variant calling</h4>
<p>Run the HaplotypeCaller on each sample's BAM file(s) (if a sample's data is spread over more than one BAM, then pass them all in together) to create single-sample gVCFs, with the option <code>--emitRefConfidence GVCF</code>, and using the <code>.g.vcf</code> extension for the output file.</p>
<p>Note that versions older than 3.4 require passing the options <code>--variant_index_type LINEAR --variant_index_parameter 128000</code> to set the correct index strategy for the output gVCF. </p>
<h4>2. Optional data aggregation step</h4>
<p>If you have more than a few hundred samples, run CombineGVCFs on batches of ~200 gVCFs to hierarchically merge them into a single gVCF. This will make the next step more tractable.</p>
<h4>3. Joint genotyping</h4>
<p>Take the outputs from step 2 (or step 1 if dealing with fewer samples) and run GenotypeGVCFs on all of them together to create the raw SNP and indel VCFs that are usually emitted by the callers.</p>
<h4>4. Variant recalibration</h4>
<p>Finally, resume the classic GATK Best Practices workflow by running VQSR on these &quot;regular&quot; VCFs according to our usual recommendations.</p>
<p>That's it! Fairly simple in practice, but we predict this is going to have a huge impact in how people perform variant discovery in large cohorts. We certainly hope it helps people deal with the challenges posed by ever-growing datasets. </p>
<p>As always, we look forward to comments and observations from the research community!</p>