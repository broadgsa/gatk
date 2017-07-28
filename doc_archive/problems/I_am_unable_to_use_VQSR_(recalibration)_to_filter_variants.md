## I am unable to use VQSR (recalibration) to filter variants

http://gatkforums.broadinstitute.org/gatk/discussion/3225/i-am-unable-to-use-vqsr-recalibration-to-filter-variants

<h3>The problem:</h3>
<p>Our preferred method for filtering variants after the calling step is to use VQSR, a.k.a. recalibration. However, it requires well-curated training/truth resources, which are typically not available for organisms other than humans, and it also requires a large amount of variant sites to operate properly, so it is not suitable for some small-scale experiments such as targeted gene panels or exome studies with fewer than 30 exomes. For the latter, it is sometimes possible to pad your cohort with exomes from another study (especially for humans -- use 1000 Genomes or ExAC!) but again for non-human organisms it is often not possible to do this. </p>
<hr />
<h3>The solution: hard-filtering</h3>
<p>So, if this is your case and you are sure that you cannot use VQSR, then you will need to use the <a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php">VariantFiltration</a> tool to <strong>hard-filter your variants</strong>. To do this, you will need to compose filter expressions using JEXL as explained <a href="http://www.broadinstitute.org/gatk/guide/article?id=1255">here</a> based on the <strong>generic filter recommendations detailed below</strong>. There is a <a href="https://www.broadinstitute.org/gatk/guide/article?id=2806">tutorial</a> that shows how to achieve this step by step. Be sure to also read the documentation explaining <a href="https://www.broadinstitute.org/gatk/guide/article?id=6925">how to understand and improve upon the generic hard filtering recommendations</a>.</p>
<hr />
<h3>But first, some caveats</h3>
<p>Let's be painfully clear about this: there is no magic formula that will give you perfect results. Filtering variants manually, using thresholds on annotation values, is subject to all sorts of caveats. The appropriateness of both the annotations and the threshold values is very highly dependent on the specific callset, how it was called, what the data was like, what organism it belongs to, etc.</p>
<p>HOWEVER, because we want to help and people always say that something is better than nothing (not necessarily true, but let's go with that for now), we have formulated some <strong>generic recommendations</strong> that should at least provide <strong>a starting point for people to experiment with their data</strong>. </p>
<p>In case you didn't catch that bit in bold there, we're saying that you absolutely SHOULD NOT expect to run these commands and be done with your analysis. You absolutely SHOULD expect to have to evaluate your results critically and TRY AGAIN with some parameter adjustments until you find the settings that are right for your data. </p>
<p>In addition, please note that these recommendations are mainly designed for dealing with very small data sets (in terms of both number of samples or size of targeted regions). If you are not using VQSR because you do not have training/truth resources available for your organism, then you should expect to have to do even more tweaking on the filtering parameters.</p>
<hr />
<h3>Filtering recommendations</h3>
<p>Here are some recommended arguments to use with VariantFiltration when ALL other options are unavailable to you. Be sure to read the documentation explaining <a href="https://www.broadinstitute.org/gatk/guide/article?id=6925">how to understand and improve upon these recommendations</a>. </p>
<p>Note that these JEXL expressions will tag as filtered any sites where the annotation value <strong>matches</strong> the expression. So if you use the expression <code>QD &lt; 2.0</code>, any site with a QD lower than 2 will be tagged as failing that filter. </p>
<h4>For SNPs:</h4>
<ul>
<li><code>QD &lt; 2.0</code></li>
<li><code>MQ &lt; 40.0</code></li>
<li><code>FS &gt; 60.0</code></li>
<li><code>SOR &gt; 3.0</code>  </li>
<li><code>MQRankSum &lt; -12.5</code></li>
<li><code>ReadPosRankSum &lt; -8.0</code></li>
</ul>
<p>If your callset was generated with UnifiedGenotyper for legacy reasons, you can add <code>HaplotypeScore &gt; 13.0</code>. </p>
<h4>For indels:</h4>
<ul>
<li><code>QD &lt; 2.0</code></li>
<li><code>ReadPosRankSum &lt; -20.0</code></li>
<li><code>InbreedingCoeff &lt; -0.8</code></li>
<li><code>FS &gt; 200.0</code></li>
<li><code>SOR &gt; 10.0</code>  </li>
</ul>
<hr />
<h3>And now some more IMPORTANT caveats (don't skip this!)</h3>
<ul>
<li>
<p>The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.</p>
</li>
<li>
<p>For shallow-coverage (&lt;10x), it is virtually impossible to use manual filtering to reliably separate true positives from false positives. You really, really, really should use the protocol involving variant quality score recalibration. If you can't do that, maybe you need to take a long hard look at your experimental design. In any case you're probably in for a world of pain.</p>
</li>
<li>The maximum DP (depth) filter only applies to whole genome data, where the probability of a site having exactly N reads given an average coverage of M is a well-behaved function.  First principles suggest this should be a binomial sampling but in practice it is more a Gaussian distribution.  Regardless, the DP threshold should be set a 5 or 6 sigma from the mean coverage across all samples, so that the DP &gt; X threshold eliminates sites with excessive coverage caused by alignment artifacts.  Note that <strong>for exomes, a straight DP filter shouldn't be used</strong> because the relationship between misalignments and depth isn't clear for capture data. </li>
</ul>
<hr />
<h3>Finally, a note of hope</h3>
<p>Some bits of this article may seem harsh, or depressing. Sorry. We believe in giving you the cold hard truth. </p>
<p>HOWEVER, we do understand that this is one of the major points of pain that GATK users encounter -- along with understanding how VQSR works, so really, whichever option you go with, you're going to suffer. </p>
<p>And we do genuinely want to help. So although we can't look at every single person's callset and give an opinion on how it looks (no, seriously, don't ask us to do that), we do want to hear from you about how we can best help you help yourself. What information do you feel would help you make informed decisions about how to set parameters? Are the meanings of the annotations not clear? Would knowing more about how they are computed help you understand how you can use them? Do you want more math? Less math, more concrete examples? </p>
<p>Tell us what you'd like to see here, and we'll do our best to make it happen. (no unicorns though, we're out of stock)</p>
<p>We also welcome testimonials from you. We are one small team; you are a legion of analysts all trying different things. Please feel free to come forward and share your findings on what works particularly well in your hands. </p>