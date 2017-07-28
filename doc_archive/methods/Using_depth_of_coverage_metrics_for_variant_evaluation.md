## Using depth of coverage metrics for variant evaluation

http://gatkforums.broadinstitute.org/gatk/discussion/4721/using-depth-of-coverage-metrics-for-variant-evaluation

<h3>Overview</h3>
<p>This document describes the proper use of metrics associated with depth of coverage for the purpose of evaluating variants.</p>
<p>The metrics involved are the following:</p>
<ul>
<li><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample (AD)</a>:</strong> outputs the depth of coverage of each allele per sample.  </li>
<li><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage (DP)</a>:</strong> outputs the filtered depth of coverage for each sample and the unfiltered depth of coverage across all samples.</li>
</ul>
<p>For an overview of the tools and concepts involved in performing sequence coverage analysis, where the purpose is to answer the common question: &quot;(Where) Do I have enough sequence data to be empowered to discover variants with reasonable confidence?&quot;, please see <a href="https://www.broadinstitute.org/gatk/guide/article?id=40">this document</a>.</p>
<hr />
<h3>Coverage annotations: DP and AD</h3>
<p>The variant callers generate two main coverage annotation metrics: the allele depth per sample (AD) and overall depth of coverage (DP, available both per sample and across all samples, with important differences), controlled by the following annotator modules:</p>
<ul>
<li><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample (AD)</a>:</strong> outputs the depth of coverage of each allele per sample.  </li>
<li><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage (DP)</a>:</strong> outputs the filtered depth of coverage for each sample and the unfiltered depth of coverage across all samples.</li>
</ul>
<p>At the sample level, these annotations are highly complementary metrics that provide two important ways of thinking about the depth of the data available for a given sample at a given site. The key difference is that the AD metric is based on unfiltered read counts while the sample-level DP is based on filtered read counts (see tool documentation for a list of read filters that are applied by default for each tool). As a result, they should be interpreted differently. </p>
<p>The sample-level DP is in some sense reflective of the power I have to determine the genotype of the sample at this site, while the AD tells me how many times I saw each of the REF and ALT alleles in the reads, free of any bias potentially introduced by filtering the reads. If, for example, I believe there really is a an A/T polymorphism at a site, then I would like to know the counts of A and T bases in this sample, even for reads with poor mapping quality that would normally be excluded from the statistical calculations going into GQ and QUAL.</p>
<p>Note that because the AD includes reads and bases that were filtered by the caller (and in case of indels, is based on a statistical computation), it  should not be used to make assumptions about the genotype that it is associated with. Ultimately, the phred-scaled genotype likelihoods (PLs) are what determines the genotype calls.</p>
<hr />
<p>TO BE CONTINUED... </p>