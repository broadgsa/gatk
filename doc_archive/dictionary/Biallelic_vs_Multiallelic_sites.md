## Biallelic vs Multiallelic sites

http://gatkforums.broadinstitute.org/gatk/discussion/6455/biallelic-vs-multiallelic-sites

<p>A <strong>biallelic</strong> site is a specific locus in a genome that contains two observed alleles, counting the reference as one, and therefore allowing for one variant allele. In practical terms, this is what you would call a site where, across multiple samples in a cohort, you have evidence for a single non-reference allele. Shown below is a toy example in which the consensus sequence for samples 1-3 have a <em>deletion</em> at position 7. Sample 4 matches the reference. This is considered a biallelic site because there are only two possible alleles-- a deletion, or the reference allele <code>G</code>.</p>
<pre><code>           1 2 3 4 5 6 7 8 9
Reference: A T A T A T G C G
Sample 1 : A T A T A T - C G
Sample 2 : A T A T A T - C G
Sample 3 : A T A T A T - C G
Sample 4 : A T A T A T G C G</code></pre>
<hr />
<p>A <strong>multiallelic</strong> site is a specific locus in a genome that contains three or more observed alleles, again counting the reference as one, and therefore allowing for two or more variant alleles. This is what you would call a site where, across multiple samples in a cohort, you see evidence for two or more non-reference alleles. Show below is a toy example in which the consensus sequences for samples 1-3 have a <em>deletion</em> or a <em>SNP</em> at the 7th position. Sample 4 matches the reference. This is considered a multiallelic site because there are four possible alleles-- a deletion, the reference allele <code>G</code>, a <code>C</code> (SNP), or a <code>T</code> (SNP). True multiallelic sites are not observed very frequently unless you look at very large cohorts, so they are often taken as a sign of a noisy region where artifacts are likely. </p>
<pre><code>           1 2 3 4 5 6 7 8 9
Reference: A T A T A T G C G
Sample 1 : A T A T A T - C G
Sample 2 : A T A T A T C C G
Sample 3 : A T A T A T T C G
Sample 4 : A T A T A T G C G</code></pre>