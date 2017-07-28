## What types of variants can GATK tools detect / handle?

http://gatkforums.broadinstitute.org/gatk/discussion/3682/what-types-of-variants-can-gatk-tools-detect-handle

<p>The answer depends on what tool we're talking about, and whether we're considering variant discovery or variant manipulation.</p>
<h4>Variant manipulation</h4>
<p>GATK variant manipulation tools are able to recognize the following types of alleles:</p>
<ul>
<li>SNP (single nucleotide polymorphism)</li>
<li>INDEL (insertion/deletion)</li>
<li>MIXED (combination of SNPs and indels at a single position)</li>
<li>MNP (multi-nucleotide polymorphism, e.g. a dinucleotide substitution)</li>
<li>SYMBOLIC (such as the <code>&lt;NON-REF&gt;</code> allele used in GVCFs produced by HaplotypeCaller, the <code>*</code> allele used to signify the presence of a <a href="https://www.broadinstitute.org/gatk/guide/article?id=6926">spanning deletion</a>, or undefined events like a very large allele or one that's fuzzy and not fully modeled; i.e. there's some event going on here but we don't know what exactly)</li>
</ul>
<p>Note that SelectVariants, the GATK tool most used for VCF subsetting operations, discriminates strictly between these categories. This means that if you use for example <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php#--selectTypeToInclude"><code>-selectType</code></a> <code>INDEL</code> to pull out indels, it will only select pure INDEL records, excluding any MIXED records that might include a SNP allele in addition to the insertion or deletion alleles of interest. To include those you would have to also specify <code>selectType MIXED</code> in the same command. </p>
<h4>Variant discovery</h4>
<p>The HaplotypeCaller is a sophisticated variant caller that can call different types of variants at the same time. So in addition to SNPs and indels, it is capable of emitting mixed records by default, as well as symbolic representations for e.g. spanning deletions. It does emit physical phasing information, but in its current version, HC is not able to emit MNPs. If you would like to combine contiguous SNPs into MNPs, you will need to use the ReadBackedPhasing tool with the MNP merging function activated. See the tool documentation for details. Our older (and now deprecated) variant caller, UnifiedGenotyper, was even more limited. It only called SNPs and indels, and did so separately (even if you ran in calling mode BOTH, the program performed separate calling operations internally) so it was not able to recognize that SNPs and Indels should be emitted together as a joint record when they occur at the same site.</p>
<p>The general release version of GATK is currently not able to detect SVs (structural variations) or CNVs (copy number variations). However, the alpha version of GATK 4 (the next generation of GATK tools) includes tools for performing CNV (copy number variation) analysis in exome data. Let us know if you're interested in trying them out by commenting on this article in the forum.</p>
<p>There is also a third-party software package called <a href="http://www.broadinstitute.org/gatk/guide/topic?name=third-party-tools">GenomeSTRiP</a> built on top of GATK that provides SV (structural variation) analysis capabilities.</p>