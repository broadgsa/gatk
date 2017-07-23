## Genotype Refinement workflow

http://gatkforums.broadinstitute.org/gatk/discussion/4723/genotype-refinement-workflow

<h3>Overview</h3>
<p>This document describes the purpose and general principles of the Genotype Refinement workflow. For the mathematical details of the methods involved, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4726">Genotype Refinement math</a> documentation. For step-by-step instructions on how to apply this workflow to your data, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4727">Genotype Refinement tutorial</a>.</p>
<hr />
<h2>1. Introduction</h2>
<p>The core GATK Best Practices workflow has historically focused on variant discovery --that is, the existence of genomic variants in one or more samples in a cohorts-- and consistently delivers high quality results when applied appropriately. However, we know that the quality of the individual genotype calls coming out of the variant callers can vary widely based on the quality of the BAM data for each sample.  The goal of the Genotype Refinement workflow is to use additional data to improve the accuracy of genotype calls and to filter genotype calls that are not reliable enough for downstream analysis. In this sense it serves as an optional extension of the variant calling workflow, intended for researchers whose work requires high-quality identification of individual genotypes.</p>
<p><strong>A few commonly asked questions are:</strong></p>
<h3>What studies can benefit from the Genotype Refinement workflow?</h3>
<p>While every study can benefit from increased data accuracy, this workflow is especially useful for analyses that are concerned with how many copies of each variant an individual has (e.g. in the case of loss of function) or with the transmission (or de novo origin) of a variant in a family.</p>
<h3>What additional data do I need to run the Genotype Refinement workflow?</h3>
<p>If  a “gold standard” dataset for SNPs is available, that can be used as a very powerful set of priors on the genotype likelihoods in your data. For analyses involving families, a pedigree file describing the relatedness of the trios in your study will provide another source of supplemental information. If neither of these applies to your data, the samples in the dataset itself can provide some degree of genotype refinement (see section 5 below for details).</p>
<h3>Is the Genotype Refinement workflow going to change my data? Can I still use my old analysis pipeline?</h3>
<p>After running the Genotype Refinement workflow, several new annotations will be added to the INFO and FORMAT fields of your variants (see below), GQ fields will be updated, and genotype calls may be modified. However, the Phred-scaled genotype likelihoods (PLs) which indicate the original genotype call (the genotype candidate with PL=0) will remain untouched. Any analysis that made use of the PLs will produce the same results as before.</p>
<hr />
<h2>2. The Genotype Refinement workflow</h2>
<h3>Overview</h3>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/40/a04b64d7884a07d1b562ead4002c60.jpg" />
<h3>Input</h3>
<p>Begin with recalibrated variants from VQSR at the end of the best practices pipeline. The filters applied by VQSR will be carried through the Genotype Refinement workflow.</p>
<h3>Step 1: Derive posterior probabilities of genotypes</h3>
<h4>Tool used: CalculateGenotypePosteriors</h4>
<p>Using the Phred-scaled genotype likelihoods (PLs) for each sample, prior probabilities for a sample taking on a HomRef, Het, or HomVar genotype are applied to derive the posterior probabilities of the sample taking on each of those genotypes. A sample’s PLs were calculated by HaplotypeCaller using only the reads for that sample. By introducing additional data like the allele counts from the 1000 Genomes project and the PLs for other individuals in the sample’s pedigree trio, those estimates of genotype likelihood can be improved based on what is known about the variation of other individuals.</p>
<p>SNP calls from the 1000 Genomes project capture the vast majority of variation across most human populations and can provide very strong priors in many cases. At sites where most of the 1000 Genomes samples are homozygous variant with respect to the reference genome, the probability of a sample being analyzed of also being homozygous variant is very high.</p>
<p>For a sample for which both parent genotypes are available, the child’s genotype can be supported or invalidated by the parents’ genotypes based on Mendel’s laws of allele transmission. Even the confidence of the parents’ genotypes can be recalibrated, such as in cases where the genotypes output by HaplotypeCaller are apparent Mendelian violations.</p>
<h3>Step 2: Filter low quality genotypes</h3>
<h4>Tool used: VariantFiltration</h4>
<p>After the posterior probabilities are calculated for each sample at each variant site, genotypes with GQ &lt; 20 based on the posteriors are filtered out. GQ20 is widely accepted as a good threshold for genotype accuracy, indicating that there is a 99% chance that the genotype in question is correct. Tagging those low quality genotypes indicates to researchers that these genotypes may not be suitable for downstream analysis.  However, as with the VQSR, a filter tag is applied, but the data is not removed from the VCF.</p>
<h3>Step 3: Annotate possible de novo mutations</h3>
<h4>Tool used: VariantAnnotator</h4>
<p>Using the posterior genotype probabilities, possible de novo mutations are tagged. Low confidence de novos have child GQ &gt;= 10 and AC &lt; 4 or AF &lt; 0.1%, whichever is more stringent for the number of samples in the dataset. High confidence de novo sites have all trio sample GQs &gt;= 20 with the same AC/AF criterion.</p>
<h3>Step 4: Functional annotation of possible biological effects</h3>
<h4>Tool options: SnpEff or Oncotator (both are non-GATK tools)</h4>
<p>Especially in the case of de novo mutation detection, analysis can benefit from the functional annotation of variants to restrict variants to exons and surrounding regulatory regions. The GATK currently does not feature integration with any functional annotation tool, but SnpEff and Oncotator are useful utilities that can work with the GATK's VCF output.</p>
<hr />
<h2>3. Output annotations</h2>
<p>The Genotype Refinement Pipeline adds several new info- and format-level annotations to each variant. GQ fields will be updated, and genotypes calculated to be highly likely to be incorrect will be changed. The Phred-scaled genotype likelihoods (PLs) carry through the pipeline without being changed. In this way, PLs can be used to derive the original genotypes in cases where sample genotypes were changed.</p>
<h3>Population Priors</h3>
<p>New INFO field annotation PG is a vector of the Phred-scaled prior probabilities of a sample at that site being HomRef, Het, and HomVar. These priors are based on the input samples themselves along with data from the supporting samples if the variant in question overlaps another in the supporting dataset.</p>
<h3>Phred-Scaled Posterior Probability</h3>
<p>New FORMAT field annotation PP is the Phred-scaled posterior probability of the sample taking on each genotype for the given variant context alleles. The PPs represent a better calibrated estimate of genotype probabilities than the PLs are recommended for use in further analyses instead of the PLs.</p>
<h3>Genotype Quality</h3>
<p>Current FORMAT field annotation GQ is updated based on the PPs. The calculation is the same as for GQ based on PLs.</p>
<h3>Joint Trio Likelihood</h3>
<p>New FORMAT field annotation JL is the Phred-scaled joint likelihood of the posterior genotypes for the trio being incorrect. This calculation is based on the PLs produced by HaplotypeCaller (before application of priors), but the genotypes used come from the posteriors. The goal of this annotation is to be used in combination with JP to evaluate the improvement in the overall confidence in the trio’s genotypes after applying CalculateGenotypePosteriors. The calculation of the joint likelihood is given as:</p>
<p>$$ -10<em>\log ( 1-GL_{mother}[\text{Posterior mother GT}] </em> GL<em>{father}[\text{Posterior father GT}] * GL</em>{child}[\text{Posterior child GT}] ) $$</p>
<p>where the GLs are the genotype likelihoods in [0, 1] probability space.</p>
<h3>Joint Trio Posterior</h3>
<p>New FORMAT field annotation JP is the Phred-scaled posterior probability of the output posterior genotypes for the three samples being incorrect. The calculation of the joint posterior is given as: </p>
<p>$$ -10<em>\log (1-GP_{mother}[\text{Posterior mother GT}] </em> GP<em>{father}[\text{Posterior father GT}] * GP</em>{child}[\text{Posterior child GT}] )$$</p>
<p>where the GPs are the genotype posteriors in [0, 1] probability space.</p>
<h3>Low Genotype Quality</h3>
<p>New FORMAT field filter lowGQ indicates samples with posterior GQ less than 20. Filtered samples tagged with lowGQ are not recommended for use in downstream analyses.</p>
<h3>High and Low Confidence De Novo</h3>
<p>New INFO field annotation for sites at which at least one family has a possible de novo mutation. Following the annotation tag is a list of the children with de novo mutations. High and low confidence are output separately.</p>
<hr />
<h2>4. Example</h2>
<p>Before:</p>
<pre><code class="pre_md">1       1226231 rs13306638      G       A       167563.16       PASS    AC=2;AF=0.333;AN=6;…        GT:AD:DP:GQ:PL  0/0:11,0:11:0:0,0,249   0/0:10,0:10:24:0,24,360 1/1:0,18:18:60:889,60,0</code class="pre_md"></pre>
<p>After:</p>
<pre><code class="pre_md">1       1226231 rs13306638      G       A       167563.16       PASS    AC=3;AF=0.500;AN=6;…PG=0,8,22;…    GT:AD:DP:GQ:JL:JP:PL:PP 0/1:11,0:11:49:2:24:0,0,249:49,0,287    0/0:10,0:10:32:2:24:0,24,360:0,32,439   1/1:0,18:18:43:2:24:889,60,0:867,43,0</code class="pre_md"></pre>
<p>The original call for the child (first sample) was HomRef with GQ0.  However, given that, with high confidence, one parent is HomRef and one is HomVar, we expect the child to be heterozygous at this site.  After family priors are applied, the child’s genotype is corrected and its GQ is increased from 0 to 49.  Based on the allele frequency from 1000 Genomes for this site, the somewhat weaker population priors favor a HomRef call (PG=0,8,22). The combined effect of family and population priors still favors a Het call for the child.</p>
<p>The joint likelihood for this trio at this site is two, indicating that the genotype for one of the samples may have been changed.  Specifically a low JL indicates that posterior genotype for at least one of the samples was not the most likely as predicted by the PLs. The joint posterior value for the trio is 24, which indicates that the GQ values based on the posteriors for all of the samples are at least 24. (See above for a more complete description of JL and JP.)</p>
<hr />
<h2>5. More information about priors</h2>
<p>The Genotype Refinement Pipeline uses Bayes’s Rule to combine independent data with the genotype likelihoods derived from HaplotypeCaller, producing more accurate and confident genotype posterior probabilities. Different sites will have different combinations of priors applied based on the overlap of each site with external, supporting SNP calls and on the availability of genotype calls for the samples in each trio.</p>
<h3>Input-derived Population Priors</h3>
<p>If the input VCF contains at least 10 samples, then population priors will be calculated based on the discovered allele count for every called variant.</p>
<h3>Supporting Population Priors</h3>
<p>Priors derived from supporting SNP calls can only be applied at sites where the supporting calls overlap with called variants in the input VCF. The values of these priors vary based on the called reference and alternate allele counts in the supporting VCF. Higher allele counts (for ref or alt) yield stronger priors.</p>
<h3>Family Priors</h3>
<p>The strongest family priors occur at sites where the called trio genotype configuration is a Mendelian violation. In such a case, each Mendelian violation configuration is penalized by a de novo mutation probability (currently 10-6). Confidence also propagates through a trio. For example, two GQ60 HomRef parents can substantially boost a low GQ HomRef child and a GQ60 HomRef child and parent can improve the GQ of the second parent. Application of family priors requires the child to be called at the site in question. If one parent has a no-call genotype, priors can still be applied, but the potential for confidence improvement is not as great as in the 3-sample case.</p>
<h3>Caveats</h3>
<p>Right now family priors can only be applied to biallelic variants and population priors can only be applied to SNPs. Family priors only work for trios.</p>