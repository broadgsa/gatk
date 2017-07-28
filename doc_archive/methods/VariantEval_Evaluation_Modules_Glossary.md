## VariantEval Evaluation Modules Glossary

http://gatkforums.broadinstitute.org/gatk/discussion/6309/varianteval-evaluation-modules-glossary

<h3>Table of Contents</h3>
<h4>Default modules:</h4>
<ul>
<li><strong><a href="#compoverlap">CompOverlap</a></strong>: gives concordance metrics based on the overlap between the evaluation and comparison file</li>
<li><strong><a href="#countvariants">CountVariants</a></strong>:  counts different types (SNP, insertion, complex, etc.) of variants present within your evaluation file and gives related metrics</li>
<li><strong>IndelLengthHistogram</strong>: gives a table of values for plotting a histogram of indel lengths found in your evaluated variants.</li>
<li><strong><a href="#indelsummary">IndelSummary</a></strong>: gives metrics related to insertions and deletions (count, multiallelic sites, het-hom ratios, etc.)</li>
<li><strong><a href="#multiallelicsummary">MultiallelicSummary</a></strong>: gives metrics relevant to multiallelic variant sites, including amount, ratio, and TiTv</li>
<li><strong><a href="#titvvariantevaluator">TiTvVariantEvaluator</a></strong>: gives the number and ratio of transition and transversion variants for your evaluation file, comparison file, and ancestral alleles</li>
<li><strong>ValidationReport</strong>: details the sensitivity and specificity of your callset, given follow-up validation assay data</li>
<li><strong>VariantSummary</strong>: gives a summary of metrics related to SNPs and indels
<h4>Other available modules:</h4></li>
<li><strong>MendelianViolationEvaluator</strong>: detects and counts Mendelian violations, given data from parent samples.</li>
<li><strong>PrintMissingComp</strong>: returns the number of variant sites present in your callset that were not found in the truth set.</li>
<li><strong>ThetaVariantEvaluator</strong>: computes different estimates of theta based on variant sites and genotypes</li>
<li><strong>MetricsCollection</strong>: includes all minimum metrics discussed in [this article]() (link to follow; document in progress). Runs by default if <em>CompOverlap</em>, <em>IndelSummary</em>, <em>TiTvVariantEvaluator</em>, <em>CountVariants</em>, &amp; <em>MultiallelicSummary</em> are run as well. (included in the nightly build for immediate use or in the next release of GATK)
<sub>* At the time of writing, the listed modules were present. To check modules present in your specific GATK version, use the <code>-list</code> <a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php#--evalModule">command</a>. </sub></li>
</ul>
<hr />
<h3>General</h3>
<p>Each table has a few columns of data that will be the same across multiple evaluation modules. To avoid listing them multiple times, they will be specified here</p>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/50/bab565feb85403f4b346b975f30eea.png" />
<ul>
<li><em>CompOverlap</em>- In the above example, we see the first column is the CompOverlap. This first column will always be the name of the evaluation module you are currently viewing. <em>IndelSummary</em> will say &quot;IndelSummary&quot;, <em>CountVariants</em> will say &quot;CountVariants&quot; and so on.</li>
<li><em>CompRod</em>- shows which file is being compared to the <code>eval</code> for that row.
By default, this is dbsnp, but you can specify additional comparison files using <code>-comp</code>, and name them using <code>:</code>. E.g. <code>-comp:name \path\to\file.vcf</code> where <code>name</code> is the name you wish to specify for the CompRod column and <code>\path\to\file.vcf</code> is your comparison file. If left unnamed, these additional comparison files will default to &quot;comp&quot; in the CompRod column.</li>
<li><em>EvalRod</em>- shows which file is being evaluated.
This is useful when specifying multiple <code>eval</code> files. They can be named using the <code>:</code> notation as above. When unnamed, they will default to &quot;eval&quot; in the EvalRod column.</li>
<li><em>JexlExpression</em>- a Jexl query that was applied to the file. For details on Jexl expressions, please read about them <a href="http://gatkforums.broadinstitute.org/discussion/1255/variant-filtering-methods-involving-jexl-queries">here</a></li>
<li><em>Novelty</em>- has three possible values; all, known, and novel. &quot;Novel&quot; includes anything seen exclusively in the eval that is not seen in the comp. &quot;Known&quot; includes anything seen in both the eval and the comp. &quot;All&quot; is the sum of &quot;Novel&quot; and &quot;Known&quot;.
By default, the comp used to determine novelty is dbsnp. To change this, you must specify <code>-knownName</code> with the new comparison file you have passed in.</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples</p>
<hr />
<p><a name="compoverlap"></a></p>
<h3>CompOverlap</h3>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a4/8b0d8da73e22cc4bd55c50b5c3bd40.png" />
<ul>
<li><em>nEvalVariants</em>- the number of variants in the <code>eval</code> file</li>
<li><em>novelSites</em>- the number of variants in the <code>eval</code> considered to be novel in comparison to dbsnp (same as novel row of nEvalVariants column)</li>
<li><em>nVariantsAtComp</em>- the number of variants present in <code>eval</code> that match the location of a variant in the comparison file (same as known row of nEvalVariants)</li>
<li><em>compRate</em>- nVariantsAtComp divided by nEvalVariants</li>
<li><em>nConcordant</em>- the number of variants present in <code>eval</code> that exactly match the genotype present in the comparison file</li>
<li><em>concordantRate</em>- nConcordant divided by nVariantsAtComp</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples </p>
<hr />
<p><a name="countvariants"></a></p>
<h3>CountVariants</h3>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/09/3d8f13f0f75ee7a67194fbd2100c4a.png" />
<ul>
<li><em>nProcessedLoci</em>- the number of loci iterated over in the reference file (also found in <em>MultiallelicSummary</em>)</li>
<li><em>nCalledLoci</em>- the number of loci called in the <code>eval</code> file</li>
<li><em>nRefLoci</em>- the number of loci in <code>eval</code> that matched the reference file</li>
<li><em>nVariantLoci</em>- the number of loci in <code>eval</code> that did not match the reference file</li>
<li><em>variantRate</em>- nVariantLoci divided by nProcessedLoci</li>
<li><em>variantRatePerBp</em>- nProcessedLoci divided by nVariantLoci (a truncated integer)</li>
<li><em>nSNPs</em>- the number of variants determined to be single-nucleotide polymorphisms</li>
<li><em>nMNPs</em>- the number of variants determined to be multi-nucleotide polymorphisms</li>
<li><em>nInsertions</em>- the number of variants determined to be insertions</li>
<li><em>nDeletions</em>- the number of variants determined to be deletions</li>
<li><em>nComplex</em>- the number of variants determined to be complex (both insertions and deletions)</li>
<li><em>nSymbolic</em>- the number of variants determined to be symbolic </li>
<li><em>nMixed</em>- the number of variants determined to be mixed (cannot be determined to be SNPs, MNPs, or indels)</li>
<li><em>nNoCalls</em>- the number of sites at which there was no variant call made</li>
<li><em>nHets</em>- the number of heterozygous loci</li>
<li><em>nHomRef</em>- the number of homozygous reference loci</li>
<li><em>nHomVar</em>- the number of homozygous variant loci</li>
<li><em>nSingletons</em>- the number of variants determined to be singletons (occur only once)</li>
<li><em>nHomDerived</em>- the number of homozygous derived variants; an ancestor had a variant at that site, but the descendant in question no longer has a variant at that site and is now homozygous reference.</li>
<li><em>heterozygosity</em>- nHets divided by nProcessedLoci</li>
<li><em>heterozygosityPerBp</em>- nProcessedLoci divided by nHets (a truncated integer)</li>
<li><em>hetHomRatio</em>- nHets divided by nHomVar</li>
<li><em>indelRate</em>- nInsertions plus nDeletions plus nComplex all divided by nProcessedLoci</li>
<li><em>indelRatePerBp</em>- nProcessedLoci divided by the sum of nInsertions, nDeletions, and nComplex (a truncated integer)</li>
<li><em>insertionDeletionRatio</em>- nInsertions divided by nDeletions</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples </p>
<hr />
<p><a name="indelsummary"></a></p>
<h3>IndelSummary</h3>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/88/6be6ad9c2dceb123741c70f4f406f1.png" />
<ul>
<li><em>n_SNPs</em>- the number of SNPs (multiallelic SNPs are counted once for each allele)</li>
<li><em>n_singleton_SNPs</em>- the number of SNP singleton loci (SNPs seen only once)</li>
<li><em>n_indels</em>- the number of indels (multiallelic indels are counted once for each allele)</li>
<li><em>n_singleton_indels</em>- the number of indel singleton loci (indels seen only once)</li>
<li><em>n_indels_matching_gold_standard</em>- the number of indel loci that match indels in the gold standard (must pass in a <code>-gold</code> parameter)</li>
<li><em>gold_standard_matching_rate</em>- n_indels_matching_gold_standard divided by n_indels</li>
<li><em>n_multiallelic_indel_sites</em>- the number of indel sites that are multiallelic</li>
<li><em>percent_of_sites_with_more_than_2_alleles</em>- n_multiallelic_indel_sites divided by the total number of indel sites</li>
<li><em>SNP_to_indel_ratio</em>- n_SNPs divided by n_indels</li>
<li><em>SNP_to_indel_ratio_for_singletons</em>- n_singleton_SNPs divided by n_singleton_indels</li>
<li><em>n_novel_indels</em>- number of indels considered to be novel in comparison to dbsnp (the novel row of the n_indels column gives the same information)</li>
<li><em>indel_novelty_rate</em>- n_novel_indels divided by n_indels</li>
<li><em>n_insertions</em>- the number of insertion variants</li>
<li><em>n_deletions</em>- the number of deletion variants</li>
<li><em>insertion_to_deletion_ratio</em>- n_insertions divided by n_deletions</li>
<li><em>n_large_deletions</em>- number of deletions with a length greater than 10</li>
<li><em>n_large_insertions</em>- number of insertions with a length greater than 10</li>
<li><em>insertion_to_deletion_ratio_for_large_indels</em>- n_large_insertions divided by n_large_deletions</li>
<li><em>n_coding_indels_frameshifting</em>- the number of indels within the coding regions of the genome which cause a frameshift</li>
<li><em>n_coding_indels_in_frame</em>- the number of indels within the coding regions of the genome which do not cause a frameshift</li>
<li><em>frameshift_rate_for_coding_indels</em>- n_coding_indels_frameshifting divided by the sum of n_coding_indels_frameshifting and n_coding_indels_in_frame</li>
<li><em>SNP_het_to_hom_ratio</em>- the number of heterozygous SNPs divided by the number of homozygous variant SNPs</li>
<li><em>indel_het_to_hom_ratio</em>- the number of heterozygous indels divided by the number of homozygous variant indels</li>
<li><em>ratio_of_1_and_2_to_3_bp_insertions</em>- the sum of one and two base pair insertions divided by three base pair insertions</li>
<li><em>ratio_of_1_and_2_to_3_bp_deletions</em>- the sum of one and two base pair deletions divided by three base pair deletions</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples </p>
<hr />
<p><a name="titvvariantevaluator"></a></p>
<h3>TiTvVariantEvaluator</h3>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/71/ccd0b26f3b09c41a0fa886f23507b0.png" />
<ul>
<li><em>nTi</em>- number of transition variants in <code>eval</code> (A&harr;G or T&harr;C)</li>
<li><em>nTv</em>- number of transversion variants in <code>eval</code> (A&harr;T or G&harr;C)</li>
<li><em>tiTvRatio</em>- nTi divided by nTv</li>
<li><em>nTiInComp</em>- number of transition variants present in the comparison file</li>
<li><em>nTvInComp</em>- number of transversion variants present in the comparison file</li>
<li><em>TiTvRatioStandard</em>- nTiInComp divided by nTvInComp</li>
<li><em>nTiDerived</em>- number of transition variants derived from ancestral alleles</li>
<li><em>nTvDerived</em>- number of transversion variants derived from ancestral alleles</li>
<li><em>tiTvDerivedRatio</em>- nTiDerived divided by nTvDerived</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples </p>
<hr />
<p><a name="multiallelicsummary"></a></p>
<h3>MultiallelicSummary</h3>
<p><strong>Example Output</strong> *</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/93/53e55aeb37052c0cc1ba269836aeb9.png" />
<ul>
<li><em>nProcessedLoci</em>- number of loci iterated over in the reference file (also found in <em>CountVariants</em>)</li>
<li><em>nSNPs</em>- number of SNPs (multiallelic SNPs are only counted once overall)</li>
<li><em>nMultiSNPs</em>- number of multiallelic SNPs (again, only counted once per loci)</li>
<li><em>processedMultiSnpRatio</em>- nMultiSNPs divided by nProcessedLoci</li>
<li><em>variantMultiSnpRatio</em>- nMultiSNPs divided by nSNPs</li>
<li><em>nIndels</em>- number of indels (multiallelic indels are only counted once overall)</li>
<li><em>nMultiIndels</em>- number of multiallelic indels (again, only counted once per loci)</li>
<li><em>processedMultiIndelRatio</em>- nMultiIndels divided by nProcessedLoci</li>
<li><em>variantMultiIndelRatio</em>- nMultiIndels divided by nIndels</li>
<li><em>nTi</em>- number of transition variants at multiallelic sites</li>
<li><em>nTv</em>- number of transversion variants at multiallelic sites</li>
<li><em>TiTvRatio</em>- nTi divided by nTv</li>
<li><em>knownSNPsPartial</em>- the number of loci at which at least one allele in <code>eval</code> was found in the known comparison file (applies only to multiallelic sites)</li>
<li><em>knownSNPsComplete</em>- the number of loci at which all alleles in <code>eval</code> were also found in the known comparison file (applies only to multiallelic sites)</li>
<li><em>SNPNoveltyRate</em>- the sum of knownSNPsPartial and knownSNPsComplete divided by nMultiSNPs</li>
</ul>
<p>*Output from a rare variant association study with &gt;1500 whole genome sequenced samples </p>