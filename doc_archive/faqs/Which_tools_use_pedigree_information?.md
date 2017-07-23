## Which tools use pedigree information?

http://gatkforums.broadinstitute.org/gatk/discussion/37/which-tools-use-pedigree-information

<p>There are two types of GATK tools that are able to use pedigree (family structure) information:</p>
<h3>Tools that require a pedigree to operate</h3>
<p><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_phasing_PhaseByTransmission.php">PhaseByTransmission</a> and <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_CalculateGenotypePosteriors.php">CalculateGenotypePosterior</a> will not run without a properly formatted pedigree file. These tools are part of the Genotype Refinement workflow, which is documented <a href="https://www.broadinstitute.org/gatk/guide/article?id=4723">here</a>.</p>
<h3>Tools that are able to generate standard variant annotations</h3>
<p>The two variant callers (HaplotypeCaller and the deprecated UnifiedGenotyper) as well as VariantAnnotator and GenotypeGVCFs are all able to use pedigree information if you request an annotation that involves population structure (e.g. Inbreeding Coefficient). To be clear though, <strong>the pedigree information is not used during the variant calling process</strong>; it is only used during the annotation step at the end.</p>
<p>If you already have VCF files that were called without pedigree information, and you want to add pedigree-related annotations (e.g to use Variant Quality Score Recalibration (VQSR) with the InbreedingCoefficient as a feature annotation), don't panic. Just run the latest version of the VariantAnnotator to re-annotate your variants, requesting any missing annotations, and make sure you pass your PED file to the VariantAnnotator as well. If you forget to provide the pedigree file, the tool will run successfully but pedigree-related annotations may not be generated (this behavior is different in some older versions).</p>
<h3>About the PED format</h3>
<p>The PED files used as input for these tools are based on PLINK pedigree files. The general description can be found <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped">here</a>.</p>
<p>For these tools, the PED files must contain only the first 6 columns from the PLINK format PED file, and no alleles, like a FAM file in PLINK. </p>