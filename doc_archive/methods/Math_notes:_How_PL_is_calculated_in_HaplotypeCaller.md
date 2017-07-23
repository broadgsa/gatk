## Math notes: How PL is calculated in HaplotypeCaller

http://gatkforums.broadinstitute.org/gatk/discussion/5913/math-notes-how-pl-is-calculated-in-haplotypecaller

<p>PL is a sample-level annotation calculated by GATK variant callers such as HaplotypeCaller, recorded in the FORMAT/sample columns of variant records in VCF files. This annotation represents the normalized <a href="https://www.broadinstitute.org/gatk/guide/article?id=4260">Phred-scaled</a> likelihoods of the genotypes considered in the variant record for each sample, as described <a href="https://www.broadinstitute.org/gatk/guide/article?id=1268">here</a>. </p>
<p>This article clarifies how the PL values are calculated and how this relates to the value of the GQ field.</p>
<hr />
<h4>Contents</h4>
<ol>
<li>The basic math</li>
<li>Example and interpretation</li>
<li>Special case: non-reference confidence model (GVCF mode)</li>
</ol>
<hr />
<h3>1. The basic math</h3>
<p>The basic formula for calculating PL is:</p>
<p>$$ PL = -10 * \log{P(Genotype | Data)} $$</p>
<p>where <code>P(Genotype | Data)</code> is the conditional probability of the Genotype given the sequence Data that we have observed. The process by which we determine the value of <code>P(Genotype | Data)</code> is described in the genotyping section of the <a href="https://www.broadinstitute.org/gatk/guide/article?id=4442">Haplotype Caller documentation</a>. </p>
<p>Once we have that probability, we simply take the log of it and multiply it by -10 to put it into <a href="https://www.broadinstitute.org/gatk/guide/article?id=4260">Phred scale</a>. Then we normalize the values across all genotypes so that the PL value of the most likely genotype is 0, which we do simply by subtracting the value of the lowest PL from all the values.</p>
<p><em>The reason we like to work in <a href="https://www.broadinstitute.org/gatk/guide/article?id=4260">Phred scale</a> is because it makes it much easier to work with the very small numbers involved in these calculations. One thing to keep in mind of course is that Phred is a log scale, so whenever we need to do a division or multiplication operation (e.g. multiplying probabilities), in Phred scale this will be done as a subtraction or addition.</em></p>
<hr />
<h3>2. Example and interpretation</h3>
<p>Here’s a worked-out example to illustrate this process. Suppose we have a site where the reference allele is A, we observed one read that has a non-reference allele T at the position of interest, and we have in hand the conditional probabilities calculated by HaplotypeCaller based on that one read (if we had more reads, their contributions would be multiplied -- or in log space, added). </p>
<p><em>Please note that the values chosen for this example have been simplified and may not be reflective of actual probabilities calculated by Haplotype Caller.</em></p>
<pre><code class="pre_md"># Alleles
Reference: A
Read: T

# Conditional probabilities calculated by HC 
P(AA | Data) = 0.000001
P(AT | Data) = 0.000100
P(TT | Data) = 0.010000</code class="pre_md"></pre>
<h4>Calculate the raw PL values</h4>
<p>We want to determine the PLs of the genotype being 0/0, 0/1, and 1/1, respectively. So we apply the formula given earlier, which yields the following values:</p>
<table class="table table-striped">
<thead>
<tr>
<th>Genotype</th>
<th>A/A</th>
<th>A/T</th>
<th>T/T</th>
</tr>
</thead>
<tbody>
<tr>
<td>Raw PL</td>
<td>-10 * log(0.000001) = 60</td>
<td>-10 * log(0.000100) = 40</td>
<td>-10 * log(0.010000) = 20</td>
</tr>
</tbody>
</table>
<p>Our first observation here is that the genotype for which the conditional probability was the highest turns out to get the lowest PL value. This is expected because, as described in the <a href="https://www.broadinstitute.org/gatk/guide/article?id=1268">VCF FAQ</a>, the PL is the <em>likelihood</em> of the genotype, which means (rather unintuitively if you’re not a stats buff) it is the probability that the genotype is <strong>not</strong> correct. So, low values mean a genotype is more likely, and high values means it’s less likely.</p>
<h4>Normalize</h4>
<p>At this point we have one more small transformation to make before we emit the final PL values to the VCF: we are going to <strong>normalize</strong> the values so that the lowest PL value is zero, and the rest are scaled relative to that. Since we’re in log space, we do this simply by subtracting the lowest value, 20, from the others, yielding the following final PL values:</p>
<table class="table table-striped">
<thead>
<tr>
<th>Genotype</th>
<th>A/A</th>
<th>A/T</th>
<th>T/T</th>
</tr>
</thead>
<tbody>
<tr>
<td>Normalized PL</td>
<td>60 - 20 = 40</td>
<td>40 - 20 = 20</td>
<td>20 - 20 = 0</td>
</tr>
</tbody>
</table>
<p>We see that there is a direct relationship between the scaling of the PLs and the original probabilities: we had chosen probabilities that were each 100 times more or less likely than the next, and in the final PLs we see that the values are spaced out by a factor of 20, which is the Phred-scale equivalent of 100. This gives us a very convenient way to estimate how the numbers relate to each other -- and how reliable the genotype assignment is -- with just a glance at the PL field in the VCF record.</p>
<h4>Genotype quality</h4>
<p>We actually formalize this assessment of genotype quality in the <strong>GQ annotation</strong>, as described also in the <a href="https://www.broadinstitute.org/gatk/guide/article?id=1268">VCF FAQ</a>.The value of GQ is simply the difference between the second lowest PL and the lowest PL (which is always 0). So, in our example GQ = 20 - 0 = 20. Note that the value of GQ is capped at 99 for practical reasons, so even if the calculated GQ is higher, the value emitted to the VCF will be 99.</p>
<hr />
<h3>3. Special case: non-reference confidence model (GVCF mode)</h3>
<p>When you run HaplotypeCaller with <code>-ERC GVCF</code> to produce a gVCF, there is an additional calculation to determine the genotype likelihoods associated with the symbolic <code>&lt;NON-REF&gt;</code> allele (which represents the possibilities that remain once you’ve eliminated the REF allele and any ALT alleles that are being evaluated explicitly).</p>
<p>The PL values for any possible genotype that includes the <code>&lt;NON-REF&gt;</code> allele have to be calculated a little differently than what is explained above because HaplotypeCaller cannot directly determine the conditional probabilities of genotypes involving <code>&lt;NON-REF&gt;</code>. Instead, it uses base quality scores to model the genotype likelihoods. </p>