## Allele Depth (AD) is lower than expected

http://gatkforums.broadinstitute.org/gatk/discussion/6005/allele-depth-ad-is-lower-than-expected

<h3>The problem:</h3>
<p>You're trying to evaluate the support for a particular call, but the numbers in the DP (total depth) and AD (allele depth) fields aren't making any sense. For example, the sum of all the ADs doesn't match up to the DP, or even more baffling, the AD for an allele that was called is zero! </p>
<p>Many users have reported being confused by variant calls where there is apparently no evidence for the called allele. For example, sometimes a VCF may contain a variant call that looks like this:</p>
<pre><code class="pre_md">2 151214 . G A 673.77 . AN=2;DP=10;FS=0.000;MLEAF=0.500;MQ=56.57;MQ0=0;NCC=0;SOR=0.693 GT:AD:DP:GQ:PL 0/1:0,0:10:38:702,0,38</code class="pre_md"></pre>
<p>You can see in the Format field the AD values are 0 for both of the alleles. However, in the Info and FORMAT fields, the DP is 10. Because the DP in the INFO field is unfiltered and the DP in the FORMAT field is filtered, you know none of the reads were filtered out by the engine's built-in read filters. And if you look at the <a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php#--bamOutput">&quot;bamout&quot;</a>, you see 10 reads covering the position! So why is the VCF reporting an AD value of 0?</p>
<hr />
<h3>The explanation: uninformative reads</h3>
<p>This is not actually a bug -- the program is doing what we expect; this is an interpretation problem. The answer lies in <strong>uninformative reads</strong>. </p>
<p>We call a read “uninformative” when it passes the quality filters, but the likelihood of the most likely allele given the read is not significantly larger than the likelihood of the second most likely allele given the read. Specifically, the difference between the Phred scaled likelihoods must be greater than 0.2 to be considered significant. In other words, that means the most likely allele must be 60% more likely than the second most likely allele. </p>
<p>Let’s walk through an example to make this clearer. Let’s say we have 2 reads and 2 possible alleles at a site. All of the reads have passed HaplotypeCaller’s quality filters, and the likelihoods of the alleles given the reads are in the table below. </p>
<table class="table table-striped">
<thead>
<tr>
<th>Reads</th>
<th>Likelihood of A</th>
<th>Likelihood of T</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>3.8708e-7</td>
<td>3.6711e-7</td>
</tr>
<tr>
<td>2</td>
<td>4.9992e-7</td>
<td>2.8425e-7</td>
</tr>
</tbody>
</table>
<p><em>Note: Keep in mind that HaplotypeCaller marginalizes the likelihoods of the haplotypes given the reads to get the likelihoods of the alleles given the reads. The table above shows the likelihoods of the alleles given the reads. For additional details, please see the <a href="https://www.broadinstitute.org/gatk/guide/article?id=4441">HaplotypeCaller method documentation</a>.</em></p>
<p>Now, let’s convert the likelihoods into Phred-scaled likelihoods. To do this, we simply take the log of the likelihoods.</p>
<table class="table table-striped">
<thead>
<tr>
<th>Reads</th>
<th>Phred-scaled likelihood of A</th>
<th>Phred-scaled likelihood of T</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>-6.4122</td>
<td>-6.4352</td>
</tr>
<tr>
<td>2</td>
<td>-6.3011</td>
<td>-6.5463</td>
</tr>
</tbody>
</table>
<p>Now, we want to determine if read 1 is informative. To do this, we simply look at the Phred scaled likelihoods of the most likely allele and the second most likely allele. The Phred scaled likelihood of the most likely allele (A) is -6.4122.The Phred-scaled likelihood of the second most likely allele (T) is -6.4352. Taking the difference between the two likelihoods gives us 0.023. Because 0.023 is Less than 0.2, read 1 is considered uninformative. </p>
<p>To determine if read 2 is informative, we take -6.3011-(-6.5463). This gives us 0.2452, which is greater than 0.2. Read 2 is considered informative.</p>
<p>How does a difference of 0.2 mean the most likely allele is ~60% more likely than the second most likely allele? Well, because the likelihoods are Phred-scaled, 0.2 = 10^0.2 = 1.585 which is approximately 60% greater. </p>
<hr />
<h3>Conclusion</h3>
<p>So, now that we know the math behind determining which reads are informative, let’s look at how this affects the record output to the VCF. If a read is considered informative, it gets counted toward the AD and DP of the variant allele in the output record. If a read is considered uninformative, it is counted towards the DP, but not the AD. That way, the AD value reflects how many reads actually contributed support for a given allele at the site. We would not want to include uninformative reads in the AD value because we don’t have confidence in them. </p>
<p><em>Please note, however, that although an uninformative read is not reported in the AD, it is still used in calculations for genotyping. In future we may add an annotation to indicate counts of reads that were considered informative vs. uninformative. Let us know in the comments if you think that would be helpful.</em></p>
<p>In most cases, you will have enough coverage at a site to disregard small numbers of uninformative reads. Unfortunately, sometimes uninformative reads are the only reads you have at a site. In this case, we report the potential variant allele, but keep the AD values 0. The uncertainty at the site will be reflected in the QG and PL values.</p>