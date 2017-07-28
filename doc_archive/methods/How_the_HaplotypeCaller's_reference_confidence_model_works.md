## How the HaplotypeCaller's reference confidence model works

http://gatkforums.broadinstitute.org/gatk/discussion/4042/how-the-haplotypecallers-reference-confidence-model-works

<p>This document describes the reference confidence model applied by HaplotypeCaller to generate genomic VCFs (gVCFS), invoked by <code>-ERC GVCF</code> or <code>-ERC BP_RESOLUTION</code> (see <a href="http://www.broadinstitute.org/gatk/guide/article?id=4017">the FAQ on gVCFs</a> for format details).</p>
<p><em>Please note that this document may be expanded with more detailed information in the near future.</em></p>
<h3>How it works</h3>
<p>The mode works by assembling the reads to create potential haplotypes, realigning the reads to their most likely haplotypes, and then projecting these reads back onto the reference sequence via their haplotypes to compute alignments of the reads to the reference. For each position in the genome we have either an ALT call (via the standard calling mechanism) or we can estimate the chance that some (unknown) non-reference allele is segregating at this position by examining the realigned reads that span the reference base. At this base we perform two calculations:</p>
<ul>
<li>Estimate the confidence that no SNP exists at the site by contrasting all reads with the ref base vs all reads with any non-reference base.</li>
<li>Estimate the confidence that no indel of size &lt; X (determined by command line parameter) could exist at this site by calculating the number of reads that provide evidence against such an indel, and from this value estimate the chance that we would not have seen the allele confidently.</li>
</ul>
<p>Based on this, we emit the genotype likelihoods (<code>PL</code>) and compute the <code>GQ</code> (from the <code>PL</code>s) for the least confidence of these two models.</p>
<p>We use a symbolic allele pair, <code>&lt;NON_REF&gt;</code>, to indicate that the site is not homozygous reference, and because we have an ALT allele we can provide allele-specific <code>AD</code> and <code>PL</code> field values.</p>
<p>For details of the gVCF format, please see <a href="http://www.broadinstitute.org/gatk/guide/article?id=4017">the document that explains what is a gVCF</a>. </p>