## Can I apply the germline variant joint calling workflow to my RNAseq data?

http://gatkforums.broadinstitute.org/gatk/discussion/7363/can-i-apply-the-germline-variant-joint-calling-workflow-to-my-rnaseq-data

<p>We have <strong>not yet validated</strong> the joint genotyping methods (HaplotypeCaller in <code>-ERC GVCF</code> mode per-sample then GenotypeGVCFs per-cohort) on RNAseq data. Our standard recommendation is to process RNAseq samples individually as laid out in the RNAseq-specific documentation. </p>
<p>However, we know that a lot of people have been trying out the joint genotyping workflow on RNAseq data, and there do not seem to be any major technical problems. You are welcome to try it on your own data, with the caveat that we cannot guarantee correctness of results, and may not be able to help you if something goes wrong. Please be sure to examine your results carefully and critically.</p>
<p>If you do pursue this, you will need to pre-process your samples according to our RNA-specific documentation, then switch to the GVCF workflow at the HaplotypeCaller stage. For filtering, it will be up to you to determine whether the hard filtering or VQSR filtering method produce best results. We have not tested any of this so we cannot provide a recommendation. Be prepared to do a lot of analysis to validate the quality of your results. </p>
<p>Good luck!</p>