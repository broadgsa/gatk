## (howto) Evaluate a callset with VariantEval

http://gatkforums.broadinstitute.org/gatk/discussion/6211/howto-evaluate-a-callset-with-varianteval

<h2>Related Documents</h2>
<ul>
<li><a href="https://www.broadinstitute.org/gatk/guide/article?id=6308">Evaluating the quality of a variant callset</a></li>
<li><a href="https://www.broadinstitute.org/gatk/guide/article?id=6186">(howto) Evaluate a callset with CollectVariantCallingMetrics</a></li>
</ul>
<h2>Context</h2>
<p>This document will walk you through use of GATK's VariantEval tool. VariantEval allows for a lot of customizability, enabling an enhanced analysis of your callset through stratification, use of additional evaluation modules, and the ability to pass in multiple truth sets. Your callset consists of variants identified by earlier steps in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">GATK best practices pipeline</a>, and now requires additional evaluation to determine where your callset falls on the spectrum of &quot;perfectly identifies all true, biological variants&quot; to &quot;only identifies artifactual or otherwise unreal variants&quot;. When variant calling, we want the callset to maximize the correct calls, while minimizing false positive calls. While very robust methods, such as Sanger sequencing, can be used to individually sequence each potential variant, statistical analysis can be used to evaluate callsets instead, saving both time and money. These callset-based analyses are accomplished by comparing relevant metrics between your samples and a known truth set, such as dbSNP. Two tools exist to examine these metrics: VariantEval in GATK, and <a href="https://www.broadinstitute.org/gatk/guide/article?id=6186">CollectVariantCallingMetrics</a> in Picard. While the latter is currently used in the Broad Institute's production pipeline, the merits to each tool, as well as the basis for variant evaluation, are discussed <a href="https://www.broadinstitute.org/gatk/guide/article?id=6308">here</a>. </p>
<hr />
<h2>Example Analysis</h2>
<pre><code>java -jar GenomeAnalysisTK.jar \
-T VariantEval \
-R reference.b37.fasta \
-eval SampleVariants.vcf \
-D dbsnp_138.b37.excluding_sites_after_129.vcf \
-noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary \
-o SampleVariants_Evaluation.eval.grp</code></pre>
<p>This command specifies the tool (<a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php">VariantEval</a>),  input files, evaluation modules to be used, and an output file to write the results to. The output will be in the form of a <a href="https://www.broadinstitute.org/gatk/guide/article?id=1244">GATKReport</a>.</p>
<h3>Input Files</h3>
<ul>
<li><code>-eval</code>: a .vcf file containing your sample(s)' variant data you wish to evaluate. The example shown here uses a whole-genome sequenced rare variant association study performed on &gt;1500 samples. You can specify multiple files to evaluate with additional <code>-eval</code> lines.</li>
<li><code>-D</code>: a dbSNP .vcf to provide the tool a reference of known variants, which can be found in the <a href="https://www.broadinstitute.org/gatk/guide/article?id=1213">GATK bundle</a></li>
<li><code>-R</code>: a reference sequence .fasta</li>
</ul>
<h3>Evaluation Modules</h3>
<p>For our example command, we will simplify our analysis and examine results using the following minimum standard modules: <em>CompOverlap</em>, <em>IndelSummary</em>, <em>TiTvVariantEvaluator</em>, <em>CountVariants</em>, &amp; <em>MultiallelicSummary</em>. These modules will provide a reasonable assessment of variant qualities while reducing the computational burden in comparison to running the default modules. In the data we ran here, &gt;1500 whole-genome-sequenced samples, this improved the run time by 5 hours and 20 minutes compared to using the default modules, which equates to a 12% time reduction. In order to do this, all default modules are removed with <code>-noEV</code>, then the minimum standard modules are added back in. This tool uses only at variants that have passed all filtration steps to calculate metrics.</p>
<ul>
<li><strong><a href="http://gatkforums.broadinstitute.org/discussion/6309/varianteval-evaluation-modules-glossary#compoverlap">CompOverlap</a></strong>: gives concordance metrics based on the overlap between the evaluation and comparison file</li>
<li><strong><a href="http://gatkforums.broadinstitute.org/discussion/6309/varianteval-evaluation-modules-glossary#countvariants">CountVariants</a></strong>:  counts different types (SNP, insertion, complex, etc.) of variants present within your evaluation file and gives related metrics</li>
<li><strong><a href="http://gatkforums.broadinstitute.org/discussion/6309/varianteval-evaluation-modules-glossary#indelsummary">IndelSummary</a></strong>: gives metrics related to insertions and deletions (count, multiallelic sites, het-hom ratios, etc.)</li>
<li><strong><a href="http://gatkforums.broadinstitute.org/discussion/6309/varianteval-evaluation-modules-glossary#multiallelicsummary">MultiallelicSummary</a></strong>: gives metrics relevant to multiallelic variant sites, including amount, ratio, and TiTv</li>
<li><strong><a href="http://gatkforums.broadinstitute.org/discussion/6309/varianteval-evaluation-modules-glossary#titvvariantevaluator">TiTvVariantEvaluator</a></strong>: gives the number and ratio of transition and transversion variants for your evaluation file, comparison file, and ancestral alleles</li>
<li><strong>MetricsCollection</strong>: includes all minimum metrics discussed in this article (the one you are currently reading). Runs by default if <em>CompOverlap</em>, <em>IndelSummary</em>, <em>TiTvVariantEvaluator</em>, <em>CountVariants</em>, &amp; <em>MultiallelicSummary</em> are run as well. (included in the nightly build for immediate use or in the 3.5 release of GATK)</li>
</ul>
<h3>Example Output</h3>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c4/3621bbe334eb86208f85001de27aca.png" />
<p>Here we see an example of the table generated by the <em>CompOverlap</em> evaluation module. The field <code>concordantRate</code> is highlighted as it is one of the metrics we examine for quality checks. Each table generated by the sample call is in the attached files list at the end of this document, which you are free to browse at your leisure. </p>
<p>It is important to note the stratification by novelty, seen in this and all other tables for this example. The row for &quot;novel&quot; includes all variants that are found in <code>SampleVariants.vcf</code> but not found in the known variants file. By default, your known variants are in dbSNP. However, if you would like to specify a different known set of variants, you can pass in a <code>-comp</code> file, and call <code>-knownName</code> on it. (See the <a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php#--known_names">VariantEval tool documentation</a> for more information) The &quot;known&quot; row includes all variants found in <code>SampleVariants.vcf</code> and the known variants file. &quot;All&quot; totals the &quot;known&quot; and &quot;novel&quot; rows. This novelty stratification is done by default, but many other stratifications can be specified; see tool documentation for more information. </p>
<p>Compiled in the below table are all of the metrics taken from various tables that we will use to ascertain the quality of the analysis.</p>
<h3>Metrics Analysis</h3>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/57/4d2214e8171cc632e147a6f3b665a5.png" />
<ul>
<li>
<p><strong>concordantRate</strong>
Referring to percent concordance, this metric is found in the <em>CompOverlap</em> table. The concordance given here is site-only; for concordance which also checks the genotype, use GenotypeConcordance from <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#GenotypeConcordanceDetailMetrics">Picard</a> or <a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php">GATK</a> Your default truth set is dbSNP, though additional truth sets can be passed into VariantEval using the <code>-comp</code> argument.<em> In the case used here, we expect (and observe) a majority of overlap between <code>eval</code> and dbSNP. The latter contains a multitude of variants and is not highly regulated, so matching a high number of <code>eval</code> variants to it is quite likely.
<sub></em> Please note: As dbSNP is our default truth set (for comparison), and our default known (for novelty determination), you will see 0 in the novel concordantRate column. If you are interested in knowing the novel concordantRate, you must specify a truth set different from the set specified as known. </sub></p>
</li>
<li>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/aa/b5f5b2bfa190da7edf3164fe841f64.png" align= "right"/><strong>nSNPs/n_SNPs &amp; nIndels/n_indels</strong>
The number of SNPs are given in <em>CountVariants</em>, <em>MultiallelicSummary</em>, and <em>IndelSummary</em>; the number of indels are given in <em>MultiallelicSummary</em> and <em>IndelSummary</em>. Different numbers are seen in each table for the same metric due to the way in which each table calculates the metric. Take the example to the right: each of the four samples give their two major alleles and though all samples have a variant at this particular loci, all are slightly different in their calls, making this a <a href="http://gatkforums.broadinstitute.org/discussion/6455/biallelic-vs-multiallelic-sites#latest">multiallelic site</a>. <br>
<em>IndelSummary</em> counts all variants separately at a multiallelic site; It thus counts 2 SNPs (one T and one C) and 1 indel (a deletion) across all samples. <em>CountVariants</em> and <em>MultiallelicSummary</em>, on the other hand, count multiallelic sites as a single variant, while still counting indels and SNPs as separate variants. Thus, they count one indel and one SNP. If you wanted to stratify by sample, all the tables would agree on the numbers for samples 1, 2, &amp; 4, as they are <a href="http://gatkforums.broadinstitute.org/discussion/6455/biallelic-vs-multiallelic-sites#latest">biallelic sites</a>. Sample 3 is multiallelic, and <em>IndelSummary</em> would count 2 SNPs, whereas <em>CountVariants</em> and <em>MultiallleicSummary</em> would count 1 SNP. Though shown here on a very small scale, the same process occurs when analyzing a whole genome or exome of variants.<br>
Our resulting numbers (~56 million SNPs &amp; ~8-11 million indels) are for a cohort of &gt;1500 whole-genome sequencing samples. Therefore, although the numbers are quite large in comparison to the ~4.4 million average variants found in <a href="http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html">Nature's 2015 paper</a>, they are within reason for a large cohort of whole genome samples.</p>
</li>
<li>
<p><strong>Indel Ratio</strong>
The indel ratio is seen twice in our tables: as <code>insertion_to_deletion_ratio</code> under <em>IndelSummary</em>, and under <em>CountVariants</em> as <code>insertionDeletionRatio</code>. Each table gives a different ratio, due to the differences in calculating indels as discussed in the previous section. In our particular sample data set, filters were run to favor detection of more rare variants. Thus the indel ratios of the loci-based table (<em>IndelSummary</em>; 0.77 &amp; 0.69) are closer to the rare ratio than the expected normal.</p>
</li>
<li><strong>tiTvRatio</strong>
While the <em>MultiallelicSummary</em> table gives a value for the TiTv of multiallelic sites, we are more interested in the overall TiTv, given by the <em>TiTvVariantEvaluator</em>. The value seen here (2.10 - 2.19) are on the higher edge of acceptable (2.0-2.1), but are still within reason.</li>
</ul>
<hr />
<h2>Note on speed performance</h2>
<p>The purpose of running the analysis with the minimum standard evaluation modules is to minimize the time spent running VariantEval. Reducing the number of evaluation modules has some effects on the total runtime; depending on the additional specifications given (stratifications, multiple <code>-comp</code> files, etc.), running with the minimum standard evaluation modules can reduce the runtime by 10-30% in comparison to running the default evaluation modules. Further reducing the runtime can be achieved through <a href="https://www.broadinstitute.org/gatk/guide/article?id=1975">multithreading</a>, using the <code>-nt</code> argument.</p>