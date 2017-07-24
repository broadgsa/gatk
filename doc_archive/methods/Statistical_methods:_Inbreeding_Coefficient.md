## Statistical methods: Inbreeding Coefficient

http://gatkforums.broadinstitute.org/gatk/discussion/8032/statistical-methods-inbreeding-coefficient

<h2>Overview</h2>
<p>Although the name Inbreeding Coefficient suggests it is a measure of inbreeding, Inbreeding Coefficient measures the excess heterozygosity at a variant site. It can be used as a proxy for poor mapping (sites that have high Inbreeding Coefficients are typically locations in the genome where the mapping is bad and reads that are in the region mismatch the region because they belong elsewhere). At least 10 samples are required (preferably many more) in order for this annotation to be calculated properly.</p>
<h3>Theory</h3>
<p>The <a href="https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle">Wikipedia article about Hardy-Weinberg principle</a> includes some very helpful information on the theoretical underpinnings of the test, as Inbreeding Coefficient relies on the math behind the Hardy-Weinberg Principle.</p>
<h3>Use in GATK</h3>
<p>We calculate Inbreeding Coefficient as </p>
<p>$$ 1-\frac{ \text{# observed heterozygotes} }{ \text{# expected heterozygotes} } $$</p>
<p>The number of observed heterozygotes can be calculated from the data. The number of expected heterozygotes is <code>2pq</code>, where <code>p</code> is the frequency of the reference allele and <code>q</code> is the frequency of the alternate allele (AF). (Please see Hardy-Weinberg Principle link above).  </p>
<p>A value of 0 suggests the site is in Hardy-Weinberg Equilibrium. Negative values of Inbreeding Coefficient could mean there are too many heterozygotes and suggest a site with bad mapping. The other nice side effect is that one of the error modes in variant calling is for all calls to be heterozygous, which this metric captures nicely. This is why we recommend filtering out variants with negative Inbreeding Coefficients. Although positive values suggest too few heterozygotes, we do not recommend filtering out positive values because they could arise from admixture of different ethnic populations. </p>
<h4>Important note:</h4>
<p>Inbreeding Coefficient is not really robust to the assumption of being unrelated. We have found that relatedness does break down the assumptions Inbreeding Coefficient is based on. For family samples, it really depends on how many families and samples you have. For example, if you have 3 families, inbreeding coefficient is not going to work. But, if you have 10,000 samples and just a few families, it should be fine. Also, if you pass in a pedigree file (*.ped), it will use that information to calculate Inbreeding Coefficient only using the founders (i.e. individuals whose parents aren't in the callset), and as long as there are &gt;= 10 of those, the data should be pretty good.</p>
<hr />
<h2>Example: Inbreeding Coefficient</h2>
<p>In this example, let's say we are working with 100 human samples, and we are trying to calculate Inbreeding Coefficient at a site that has A for the reference allele and T for the alternate allele. </p>
<h3>Step 1: Count the number of samples that have each genotype</h3>
<p>HOM-REF A/A : 51
HET A/T : 11
HOM-VAR T/T : 38</p>
<h3>Step 2: Get all necessary information to solve equation</h3>
<p>We need to find the # observed hets and # expected hets:</p>
<p>$$ \text{number of observed hets} = 11 $$</p>
<p>from the number of observed A/T given above, and </p>
<p>$$ \text{number of expected hets} = 2pq * \text{total genotypes} $$</p>
<p>where <code>2pq</code> is the frequency of heterozygotes according to Hardy-Weinberg Equilibrium. </p>
<p>We need to multiply that frequency by the number of all genotypes in the population to get the expected number of heterozygotes.</p>
<p>So let's calculate <code>p</code>:</p>
<p>$$ p = \text{frequency of ref allele} = \frac{ \text{# ref alleles} }{ \text{total # alleles} } $$
$$ p = \frac{ 2 <em> 51 + 11 }{ 2 </em> 51 + 11 <em> 2 + 38 </em> 2} $$
$$ p = \frac{ 113 }{ 200 } = 0.565 $$</p>
<p>And now let's calculate <code>q</code>:</p>
<p>$$ q = \text{frequency of alt allele} = \frac{ \text{# alt alleles} }{ \text{total # alleles} } $$
$$ q = \frac{ 2 <em> 38 + 11 }{ 2 </em> 51 + 11 <em> 2 + 38 </em> 2 } $$
$$ q = 87/200 = 0.435 $$</p>
<p>Remember that homozygous genotypes have two copies of the allele of interest (because we're assuming a diploid organism).</p>
<p>$$ \text{number of expected hets} = 2pq <em> 100 $$
$$ = 2 </em> 0.565 <em> 0.435 </em> 100 = 49.155 $$</p>
<h3>Step 3: Plug in the Numbers</h3>
<p>$$ \text{Inbreeding Coefficient} = 1 - \frac{ \text{# observed hets} }{ \text{#expected hets} } $$
$$ \text{IC} = 1 - \frac{ 11 }{49.155} = 0.776 $$</p>
<h3>Step 4: Interpret the output</h3>
<p>Our Inbreeding Coefficient is 0.776. Because it is a positive number, we can see there are fewer than the expected number of heterozygotes according to the Hardy-Weinberg Principle. Too few heterozygotes can imply inbreeding. Depending on the cohort we are working with, this could be a sign of false positives.</p>