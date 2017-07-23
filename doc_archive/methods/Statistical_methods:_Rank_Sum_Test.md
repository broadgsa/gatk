## Statistical methods: Rank Sum Test

http://gatkforums.broadinstitute.org/gatk/discussion/8031/statistical-methods-rank-sum-test

<h2>Overview</h2>
<p>The Rank Sum Test, also known as Mann-Whitney-Wilcoxon U-test after its developers (who are variously credited in subsets and in different orders depending on the sources you read) is a statistical test that aims to determine whether there is significant difference in the values of two populations of data.</p>
<h3>Theory</h3>
<p>The <a href="https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test">Wikipedia article about the Rank Sum Test</a> includes some very helpful information on the theoretical underpinnings of the test, as well as various examples of how it can be applied.  </p>
<h3>Use in GATK</h3>
<p>This test is used by several GATK annotations, including two standard annotations that are used for variant recalibration in the Best Practices:  <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSum</a> and <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_ReadPosRankSumTest.php">ReadPosRankSum</a>. In all cases, the idea is to check, for a given candidate variant, whether the properties of the data that support the reference allele are similar to those of the data that support a variant allele. If they are not similar, we conclude that there may be some technical bias and that the candidate variant may be an artifact. </p>
<hr />
<h2>Example: BaseQualityRankSumTest</h2>
<p><em>Note: this example applies Method 2 from the Wikipedia article linked above.</em></p>
<p>In this example, we have a set of 20 reads, 10 of which support the reference allele and 10 of which support the alternate allele. At first glance, that looks like a clear heterozygous 0/1 site. But to be thorough in our analysis and to account for any technical bias, we want to determine if there is a significant difference in the base qualities of the bases that support the reference allele vs. the bases that support the alternate allele. </p>
<p>Before we proceed, we must define our null hypothesis and alternate hypothesis. </p>
<p>-<em>Null hypothesis:</em> There is <strong>no</strong> difference in the base qualities that support the reference allele and the base qualities that support the alternate allele.</p>
<p>-<em>Alternate hypothesis:</em> There <strong>is</strong> a difference in the base qualities that support the reference allele and the base qualities that support the alternate allele.</p>
<h3>Step 1: List the relevant observations</h3>
<p>Reference allele base qualities: 20, 25, 26, 30, 32, 40, 47, 50, 53, 60
Alternate allele base qualities: 0, 7, 10, 17, 20, 21, 30, 34, 40, 45</p>
<h3>Step 2: Rank the observations</h3>
<p>First, we arrange all the observations (base qualities) into a list of values ordered from lowest to highest (reference bases are in bold).</p>
<p>0, 7, 10, 17, <strong>20</strong>, 20, 21, <strong>25</strong>, <strong>26</strong>, <strong>30</strong>, 30, <strong>32</strong>, 34, <strong>40</strong>, 40, 45, <strong>47</strong>, <strong>50</strong>, <strong>53</strong>, <strong>60</strong></p>
<p>Next we determine the ranks of the values. Since there are 20 observations (the base qualities), we have 20 ranks to assign. Whenever there are ties between observations for the rank, we take the rank to be equal to the midpoint of the ranks. For example, for 20(ref) and 20(alt), we have a tie in values, so we assign each observation a rank of (5+6)/2 = 5.5.</p>
<p>The ranks from the above list are (reference ranks are in bold):</p>
<p>1, 2, 3, 4, <strong>5.5</strong>, 5.5, 7, <strong>8</strong>, <strong>9</strong>, <strong>10.5</strong>, 10.5, <strong>12</strong>, 13, <strong>14.5</strong>, 14.5, 16, <strong>17</strong>, <strong>18</strong>, <strong>19</strong>, <strong>20</strong></p>
<h3>Step 3: Add up the ranks for each group</h3>
<p>We now need to add up the ranks for the base qualities that came from the reference allele and the alternate allele.</p>
<p>$$ Rank_{ref} = 133.5 $$</p>
<p>$$ Rank_{alt} = 76.5 $$</p>
<h3>Step 4: Calculate U for each group</h3>
<p>U is a statistic that tells us the difference between the two rank totals. We can use the U statistic to calculate the z-score (explained below), which will give us our p-value.</p>
<p>Calculate U for each group (n = number of observations in each sample)</p>
<p>$$ U<em>{ref} = \frac{ n</em>{ref} <em> n<em>{alt} + n</em>{ref} </em> (n<em>{ref}+ 1) }{ 2 } - Rank</em>{ref} $$</p>
<p>$$ U<em>{alt} = \frac{ n</em>{alt} <em> n<em>{ref} + n</em>{alt} </em> (n<em>{alt} + 1) }{ 2 } - Rank</em>{alt} $$</p>
<p>$$ U_{ref} = \frac{ 10 <em> 10 + 10 </em> 11 }{ 2 } - 133.5 = 21.5 $$</p>
<p>$$ U_{alt} = \frac{ 10 <em> 10 + 10 </em> 11 }{ 2 } - 76.5 = 78.5 $$</p>
<h3>Step 5: Calculate the overall z-score</h3>
<p>Next, we need to calculate the z-score which will allow us to get the p-value. The z-score is a normalized score that allows us to compare the probability of the U score occurring in our distribution.
<a href="https://statistics.laerd.com/statistical-guides/standard-score.php">https://statistics.laerd.com/statistical-guides/standard-score.php</a></p>
<p>The equation to get the z-score is:</p>
<p>$$ z = \frac{U - mu}{u} $$ </p>
<p>Breaking this equation down:</p>
<p>$$ z = z-score $$</p>
<p>$$ U = \text{lowest of the U scores calculated in previous steps} $$</p>
<p>$$ mu = \text{mean of the U scores above} = \frac{ n<em>{ref} * n</em>{alt} }{ 2 } $$</p>
<p>$$ u = \text{standard deviation of U} = \sqrt{ \frac{n<em>{ref} * n</em>{alt} * (n<em>{ref} + n</em>{alt} + 1) }{ 12 } }  $$</p>
<p>To calculate our z:</p>
<p>$$ U = 21.5 $$</p>
<p>$$ mu = \frac{10 * 10 }{ 2 } = 50 $$</p>
<p>$$ u = \sqrt{ \frac{10 <em> 10 </em>(10 + 10 + 1) }{ 12 } } = 13.229 $$</p>
<p>So altogether we have: </p>
<p>$$ z = \frac{ 21.5 - 50 }{ 13.229 } = -2.154 $$</p>
<h3>Step 6: Calculate and interpret the p-value</h3>
<p>The p-value is the probability of obtaining a z-score at least as extreme as the one we got, assuming the null hypothesis is true. In our example, the p-value gives us the probability that there is no difference in the base qualities that support the reference allele and the base qualities that support the alternate allele. The lower the p-value, the less likely it is that there is no difference in the base qualities.</p>
<p>Going to the z-score table, or just using a <a href="http://graphpad.com/quickcalcs/pValue2/">p-value calculator</a>, we find the p-value to be 0.0312.</p>
<p>This means there is a .0312 chance that the base quality scores of the reference allele and alternate allele are the same. Assuming a p-value cutoff of 0.05, meaning there is less than 5% chance there is no difference in the two groups, and greater than or equal to 95% chance that there is a difference between the two groups, we have enough evidence to <strong>reject our null hypothesis</strong> that there is no difference in the base qualities of the reference and alternate allele. This indicates there is some bias and that the alternate allele is less well supported by the data than the allele counts suggest.</p>