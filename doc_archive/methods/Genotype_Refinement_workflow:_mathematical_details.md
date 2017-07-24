## Genotype Refinement workflow: mathematical details

http://gatkforums.broadinstitute.org/gatk/discussion/4726/genotype-refinement-workflow-mathematical-details

<h3>Overview</h3>
<p>This document describes the mathematical details of the methods involved in the Genotype Refinement workflow. For an explanation of the purpose and general principles involved in this workflow, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4723">main Genotype Refinement workflow article</a>. For step-by-step instructions on how to apply this workflow to your data, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4727">Genotype Refinement tutorial</a>.</p>
<hr />
<h2>1. Review of Bayes’s Rule</h2>
<p>HaplotypeCaller outputs the likelihoods of observing the read data given that the genotype is actually HomRef, Het, and HomVar. To convert these quantities to the probability of the genotype given the read data, we can use Bayes’s Rule. Bayes’s Rule dictates that the probability of a parameter given observed data is equal to the likelihood of the observations given the parameter multiplied by the prior probability that the parameter takes on the value of interest, normalized by the prior times likelihood for all parameter values:</p>
<p>$$ P(\theta|Obs) = \frac{P(Obs|\theta)P(\theta)}{\sum_{\theta} P(Obs|\theta)P(\theta)} $$</p>
<p>In the best practices pipeline, we interpret the genotype likelihoods as probabilities by implicitly converting the genotype likelihoods to genotype probabilities using non-informative or flat priors, for which each genotype has the same prior probability. However, in the Genotype Refinement Pipeline we use independent data such as the genotypes of the other samples in the dataset, the genotypes in a “gold standard” dataset, or the genotypes of the other samples in a family to construct more informative priors and derive better posterior probability estimates.</p>
<hr />
<h2>2. Calculation of Population Priors</h2>
<p>Given a set of samples in addition to the sample of interest (ideally non-related, but from the same ethnic population), we can derive the prior probability of the genotype of the sample of interest by modeling the sample’s alleles as two independent draws from a pool consisting of the set of all the supplemental samples’ alleles. (This follows rather naturally from the Hardy-Weinberg assumptions.) Specifically, this prior probability will take the form of a multinomial Dirichlet distribution parameterized by the allele counts of each allele in the supplemental population.  In the biallelic case the priors can be calculated as follows:</p>
<p>$$ P(GT = HomRef) = \dbinom{2}{0} \ln \frac{\Gamma(nSamples)\Gamma(RefCount + 2)}{\Gamma(nSamples + 2)\Gamma(RefCount)} $$</p>
<p>$$ P(GT = Het) = \dbinom{2}{1} \ln \frac{\Gamma(nSamples)\Gamma(RefCount + 1)\Gamma(AltCount + 1)}{\Gamma(nSamples + 2)\Gamma(RefCount)\Gamma(AltCount)} $$</p>
<p>$$ P(GT = HomVar) = \dbinom{2}{2} \ln \frac{\Gamma(nSamples)\Gamma(AltCount + 2)}{\Gamma(nSamples + 2)\Gamma(AltCount)} $$</p>
<p>where Γ is the <a href="http://en.wikipedia.org/wiki/Gamma_function">Gamma function</a>, an extension of the factorial function.</p>
<p>The prior genotype probabilities based on this distribution scale intuitively with number of samples. For example, a set of 10 samples, 9 of which are HomRef yield a prior probability of another sample being HomRef with about 90% probability whereas a set of 50 samples, 49 of which are HomRef yield a 97% probability of another sample being HomRef.</p>
<hr />
<h2>3. Calculation of Family Priors</h2>
<p>Given a genotype configuration for a given mother, father, and child trio, we set the prior probability of that genotype configuration as follows:</p>
<p>$$ P(G_M,G_F,G_C) = P(\vec{G}) \cases{ 1-10\mu-2\mu^2 &amp; no MV \cr \mu &amp; 1 MV \cr \mu^2 &amp; 2 MVs} $$</p>
<p>where the 10 configurations with a single Mendelian violation are penalized by the de novo mutation probability μ and the two configurations with two Mendelian violations by μ^2. The remaining configurations are considered valid and are assigned the remaining probability to sum to one.</p>
<p>This prior is applied to the joint genotype combination of the three samples in the trio. To find the posterior for any single sample, we marginalize over the remaining two samples as shown in the example below to find the posterior probability of the child having a HomRef genotype:</p>
<p>$$ P(G_C = HomRef | \vec{D}) = \frac{L_C(G<em>C = HomRef) \sum</em>{G_F,G_M} L_F(G_F)L_M(G<em>M)P(\vec{G})}{\sum</em>{\vec{H}}P(\vec{D}|\vec{H})P(\vec{H})} $$</p>
<p>This quantity P(Gc|D) is calculated for each genotype, then the resulting vector is Phred-scaled and output as the Phred-scaled posterior probabilities (PPs).</p>
<hr />
<h2>4. Order of the workflow</h2>
<p>Family priors are calculated and applied before population priors. The opposite ordering results in overly strong population priors because they are applied to the child and parents and then compounded when the trio likelihoods are multiplied together.</p>