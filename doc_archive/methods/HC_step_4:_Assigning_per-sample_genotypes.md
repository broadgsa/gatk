## HC step 4: Assigning per-sample genotypes

http://gatkforums.broadinstitute.org/gatk/discussion/4442/hc-step-4-assigning-per-sample-genotypes

<p>This document describes the procedure used by HaplotypeCaller to assign genotypes to individual samples based on the allele likelihoods calculated in the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4441">previous step</a>. For more context information on how this fits into the overall HaplotypeCaller method, please see the more general <a href="http://www.broadinstitute.org/gatk/guide/article?id=4148">HaplotypeCaller documentation</a>. See also the documentation on <a href="https://www.broadinstitute.org/gatk/guide/article?id=7258">the QUAL score</a> as well as <a href="https://www.broadinstitute.org/gatk/guide/article?id=5913">PL and GQ</a>.</p>
<p>Note that this describes the <strong>regular mode</strong> of HaplotypeCaller, which does not emit an estimate of reference confidence. For details on how the reference confidence model works and is applied in <code>-ERC</code> modes (<code>GVCF</code> and <code>BP_RESOLUTION</code>) please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4042">reference confidence model documentation</a>.</p>
<h3>Overview</h3>
<p>The previous step produced a table of per-read allele likelihoods for each candidate variant site under consideration. Now, all that remains to do is to evaluate those likelihoods in aggregate to determine what is the most likely genotype of the sample at each site. This is done by applying Bayes' theorem to calculate the likelihoods of each possible genotype, and selecting the most likely. This produces a genotype call as well as the calculation of various metrics that will be annotated in the output VCF if a variant call is emitted.</p>
<hr />
<h3>1. Preliminary assumptions / limitations</h3>
<h4>Quality</h4>
<p>Keep in mind that we are trying to infer the genotype of each sample given the observed sequence data, so the degree of confidence we can have in a genotype depends on both the quality and the quantity of the available data. By definition, low coverage and low quality will both lead to lower confidence calls. The GATK only uses reads that satisfy certain mapping quality thresholds, and only uses “good” bases that satisfy certain base quality thresholds (see documentation for default values). </p>
<h4>Ploidy</h4>
<p>Both the HaplotypeCaller and GenotypeGVCFs (but not UnifiedGenotyper) assume that the organism of study is diploid by default, but desired ploidy can be set using the <code>-ploidy</code> argument. The ploidy is taken into account in the mathematical development of the Bayesian calculation. The generalized form of the genotyping algorithm that can handle ploidies other than 2 is available as of version 3.3-0. Note that using ploidy for pooled experiments is subject to some practical limitations due to the number of possible combinations resulting from the interaction between ploidy and the number of alternate alleles that are considered (currently, the maximum &quot;workable&quot; ploidy is ~20 for a max number of alt alleles = 6). Future developments will aim to mitigate those limitations.</p>
<h4>Paired end reads</h4>
<p>Reads that are mates in the same pair are not handled together in the reassembly, but if they overlap, there is some special handling to ensure they are not counted as independent observations. </p>
<h4>Single-sample vs multi-sample</h4>
<p>We apply different genotyping models when genotyping a single sample as opposed to multiple samples together (as done by HaplotypeCaller on multiple inputs or GenotypeGVCFs on multiple GVCFs). The multi-sample case is not currently documented for the public but is an extension of previous work by Heng Li and others. </p>
<hr />
<h3>2. Calculating genotype likelihoods using Bayes' Theorem</h3>
<p>We use the approach described in <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/">Li 2011</a> to calculate the posterior probabilities of non-reference alleles (Methods 2.3.5 and 2.3.6) extended to handle multi-allelic variation. </p>
<p>The basic formula we use for all types of variation under consideration (SNPs, insertions and deletions) is:</p>
<p>$$ P(G|D) = \frac{ P(G) P(D|G) }{ \sum_{i} P(G_i) P(D|G_i) } $$</p>
<p>If that is meaningless to you, please don't freak out -- we're going to break it down and go through all the components one by one. First of all, the term on the left:</p>
<p>$$ P(G|D) $$</p>
<p>is the quantity we are trying to calculate for each possible genotype: the conditional probability of the genotype <strong>G</strong> given the observed data <strong>D</strong>. </p>
<p>Now let's break down the term on the right: </p>
<p>$$ \frac{ P(G) P(D|G) }{ \sum_{i} P(G_i) P(D|G_i) } $$</p>
<p>We can ignore the denominator (bottom of the fraction) because it ends up being the same for all the genotypes, and the point of calculating this likelihood is to determine the most likely genotype. The important part is the numerator (top of the fraction):</p>
<p>$$ P(G) P(D|G) $$</p>
<p>which is composed of two things: the prior probability of the genotype and the conditional probability of the data given the genotype.</p>
<p>The first one is the easiest to understand. The prior probability of the genotype <strong>G</strong>:</p>
<p>$$ P(G) $$</p>
<p>represents how probably we expect to see this genotype based on previous observations, studies of the population, and so on. By default, the GATK tools use a flat prior (always the same value) but you can input your own set of priors if you have information about the frequency of certain genotypes in the population you're studying. </p>
<p>The second one is a little trickier to understand if you're not familiar with Bayesian statistics. It is called the conditional probability of the data given the genotype, but what does that mean? Assuming that the genotype <strong>G</strong> is the true genotype, </p>
<p>$$ P(D|G) $$</p>
<p>is the probability of observing the sequence data that we have in hand. That is, how likely would we be to pull out a read with a particular sequence from an individual that has this particular genotype? We don't have that number yet, so this requires a little more calculation, using the following formula:</p>
<p>$$ P(D|G) = \prod{j} \left( \frac{P(D_j | H_1)}{2} + \frac{P(D_j | H_2)}{2} \right) $$ </p>
<p>You'll notice that this is where the diploid assumption comes into play, since here we decomposed the genotype <strong>G</strong> into:</p>
<p>$$ G = H_1H_2 $$</p>
<p>which allows for exactly two possible haplotypes. In future versions we'll have a generalized form of this that will allow for any number of haplotypes. </p>
<p>Now, back to our calculation, what's left to figure out is this:</p>
<p>$$ P(D_j|H_n) $$</p>
<p>which as it turns out is the conditional probability of the data given a particular haplotype (or specifically, a particular allele), aggregated over all supporting reads. Conveniently, that is exactly what we calculated in Step 3 of the HaplotypeCaller process, when we used the PairHMM to produce the likelihoods of each read against each haplotype, and then marginalized them to find the likelihoods of each read for each allele under consideration. So all we have to do at this point is plug the values from that table into the equation above, and we can work our way back up to obtain:</p>
<p>$$ P(G|D) $$</p>
<p>for the genotype <strong>G</strong>. </p>
<hr />
<h3>3. Selecting a genotype and emitting the call record</h3>
<p>We go through the process of calculating a likelihood for each possible genotype based on the alleles that were observed at the site, considering every possible combination of alleles. For example, if we see an A and a T at a site, the possible genotypes are AA, AT and TT, and we end up with 3 corresponding probabilities. We pick the largest one, which corresponds to the most likely genotype, and assign that to the sample. </p>
<p>Note that depending on the variant calling options specified in the command-line, we may only emit records for actual variant sites (where at least one sample has a genotype other than homozygous-reference) or we may also emit records for reference sites. The latter is discussed in the reference confidence model documentation. </p>
<p>Assuming that we have a non-ref genotype, all that remains is to calculate the various site-level and genotype-level metrics that will be emitted as annotations in the variant record, including <a href="https://www.broadinstitute.org/gatk/guide/article?id=7258">QUAL</a> as well as <a href="https://www.broadinstitute.org/gatk/guide/article?id=5913">PL and GQ</a> -- see the linked docs for details. For more information on how the other variant context metrics are calculated, please see the corresponding variant annotations documentation. </p>