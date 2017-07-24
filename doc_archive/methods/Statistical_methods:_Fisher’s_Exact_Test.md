## Statistical methods: Fisher’s Exact Test

http://gatkforums.broadinstitute.org/gatk/discussion/8056/statistical-methods-fisher-s-exact-test

<h3>Overview</h3>
<p>Fisher’s Exact Test is a statistical test that is used to analyze contingency tables, where contingency tables are matrices that contain the frequencies of the variables in play. According to statistics lore, noted statistician R.A.Fisher invented the test to determine if Dr. Muriel Bristol could actually tell the difference between milk being added to her tea or tea being added to her milk (she couldn’t). Fisher’s Exact Test is so named because it allows us to calculate the exact p-value for the experiment, rather than having to rely on an approximation. The p-value gives us the probability of observing the set of results we obtained if the null hypothesis were true, <em>i.e.</em> getting those results purely by chance. </p>
<hr />
<h3>Mathematical theory</h3>
<p>The <a href="http://mathworld.wolfram.com/FishersExactTest.html">Wolfram Math World article on Fisher’s Exact Test</a> includes some very helpful information on the theoretical underpinnings of the test, as well as an example of how it can be applied. </p>
<hr />
<h3>Use in GATK</h3>
<p>In GATK, we use Fisher’s Exact Test to calculate the FisherStrand annotation, which is an indicator of strand bias, a common source of artifactual calls. The test determines whether there is a difference in the number of reads that support the reference allele and alternate allele on each strand (<em>i.e.</em> number of reads in forward and reverse orientation). The value is reported in the FisherStrand annotation, FS in the VCF. </p>
<hr />
<h3>Example: Fisher Strand in practice</h3>
<p><em>Note: This example follows the steps given in the Wolfram article linked above.</em></p>
<p>In this example, we want to determine if there is a difference in the number of reads that support the reference allele and alternate allele on each strand. Our null hypothesis is that there is no difference in the number of reads that support the reference allele and alternate allele on each strand (there is no strand  bias). We will calculate a <a href="http://mathworld.wolfram.com/P-Value.html">p-value</a> that tells us the probability of observing our data if our null hypothesis is true (or, that there is no strand bias). The lower the p-value, the less likely we are to believe that there is no strand bias.</p>
<p>Let’s say we have 3 reads supporting the reference allele on the forward strand and 0 reads supporting the reference allele on the reverse strand. We also have 0 reads supporting the alternate allele on the forward strand and 3 reads supporting the alternate allele on the reverse strand.</p>
<p>The contingency table, or matrix, looks like this:</p>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">Forward Strand</th>
<th style="text-align: left;">Reverse Strand</th>
<th style="text-align: center;">Total</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Reference Allele</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">0</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Alternate Allele</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Total</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">6</td>
</tr>
</tbody>
</table>
<p>At first glance, it seems obvious there is some bias going on here, because each allele is only seen either on the forward strand or the reverse strand. To determine with confidence whether there really is strand bias, we will perform Fisher’s Exact Test on this set of observations.</p>
<p>We first use the <a href="http://mathworld.wolfram.com/HypergeometricDistribution.html">hypergeometric probability function</a> to calculate the probability of getting the exact matrix we have above. The probability calculation for a 2 x 2 matrix is:</p>
<p>$$P = \frac{(R<em>{1}! \times R</em>{2}! \times C<em>{1}! \times C</em>{2}!) }{ N! \times \prod<em>{ij} a</em>{ij} } $$</p>
<p>Let’s define the variables in that equation:</p>
<ul>
<li>R1 = sum of row 1</li>
<li>R2 = sum of row 2</li>
<li>C1 = sum of column 1</li>
<li>C2 = sum of column 2</li>
<li>N = R1 + R2 = C1 + C2</li>
<li>aij = values in matrix where i and j are row and column numbers</li>
</ul>
<p>Now, let’s calculate the probability P for our own matrix above:</p>
<p>$$P = \frac{3! \times 3! \times 3! \times 3!}{6! \times 3! \times 0! \times 0! \times 3!} = 0.05 $$</p>
<p>That gives us the probability of observing our own data. However, for our test, we need the probability of observing our own data <em>and</em> more extreme data. So now we need to calculate the probability of observing more extreme data, which we'll define as any matrix that has the same row and column totals as our own, and also has a probability equal to or less than our matrix probability. </p>
<h4>Matrix probability calculations</h4>
<p>Let's find all possible matrices of non-negative integers that would be consistent with the given row and column totals (i.e. total number of observations) and calculate their probability using the formula for above.</p>
<ul>
<li>Original matrix (our experimental observations)</li>
</ul>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">Forward Strand</th>
<th style="text-align: left;">Reverse Strand</th>
<th style="text-align: center;">Total</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Reference Allele</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">0</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Alternate Allele</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Total</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">6</td>
</tr>
</tbody>
</table>
<p>$$P = \frac{3! \times 3! \times 3! \times 3!}{6! \times 3! \times 0! \times 0! \times 3!} = 0.05 $$</p>
<ul>
<li>Hypothetical matrix 1</li>
</ul>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">Forward Strand</th>
<th style="text-align: left;">Reverse Strand</th>
<th style="text-align: center;">Total</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Reference Allele</td>
<td style="text-align: left;">2</td>
<td style="text-align: left;">1</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Alternate Allele</td>
<td style="text-align: left;">1</td>
<td style="text-align: left;">2</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Total</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">6</td>
</tr>
</tbody>
</table>
<p>$$P = \frac{3! \times 3! \times 3! \times 3!}{6! \times 2! \times 1! \times 1! \times 2!} = 0.45 $$</p>
<ul>
<li>Hypothetical matrix 2</li>
</ul>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">Forward Strand</th>
<th style="text-align: left;">Reverse Strand</th>
<th style="text-align: center;">Total</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Reference Allele</td>
<td style="text-align: left;">1</td>
<td style="text-align: left;">2</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Alternate Allele</td>
<td style="text-align: left;">2</td>
<td style="text-align: left;">1</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Total</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">6</td>
</tr>
</tbody>
</table>
<p>$$P = \frac{3! \times 3! \times 3! \times 3!}{6! \times 1! \times 2! \times 2! \times 1!} = 0.45 $$</p>
<ul>
<li>Hypothetical matrix 3</li>
</ul>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"></th>
<th style="text-align: left;">Forward Strand</th>
<th style="text-align: left;">Reverse Strand</th>
<th style="text-align: center;">Total</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Reference Allele</td>
<td style="text-align: left;">0</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Alternate Allele</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">0</td>
<td style="text-align: center;">3</td>
</tr>
<tr>
<td style="text-align: left;">Total</td>
<td style="text-align: left;">3</td>
<td style="text-align: left;">3</td>
<td style="text-align: center;">6</td>
</tr>
</tbody>
</table>
<p>$$P = \frac{3! \times 3! \times 3! \times 3!}{6! \times 0! \times 3! \times 3! \times 0!} = 0.05 $$</p>
<h4>Results</h4>
<p>We see that the only matrix with a probability less than or equal to our matrix is hypothetical matrix 3. We will now add the probabilities of our own matrix and matrix 3 to get the final p-value.</p>
<p>Sum all p-values less than or equal to 0.05 to calculate overall P-value:</p>
<p>$$P_{total} = 0.05\ \text{(original)} + 0.05\ \text{(matrix 3)} = 0.1 $$</p>
<p>The p-value of 0.1 tells us there is a 10% chance that there is no statistically convincing evidence of bias, despite our strong intuition that the numbers look biased. This is because there are only 6 reads, and we can’t confidently say that there is really strand bias at work based on so few reads (observations). If we had seen more, we may have had more evidence to confidently say there is bias -- or we might have realized there is no bias at this site, and the numbers we saw were an accidental effect. If you’d like to see how our confidence scales with read numbers, try working out several cases with larger numbers of reads. You’ll need to draw up a lot of possible matrices!</p>
<p>Anyway, in the GATK context we still want to transform our FS annotation value to Phred scale for convenience before writing it out to the output VCF. To get the Phred-scaled p-value, we simply plug in the p-value of 0.1 into the Phred equation like this:</p>
<p>$$ \text{Phred Score} = -10 \times \log<em>{10} \text{p-value} = -10 \times \log</em>{10} 0.1 = 10 $$</p>
<p>So the value of FS at this site would be 10. Note if we had a p-value of 1, meaning there is a 100% chance of there being no bias, the Phred score would be 0. So, a Phred-score closer to 0 means there is a lower chance of there being bias. Higher FS values therefore indicate more bias. See the documentation article on <a href="link">understanding hard-filtering recommendations</a> for more commentary on how we interpret the value of FS in practice. </p>