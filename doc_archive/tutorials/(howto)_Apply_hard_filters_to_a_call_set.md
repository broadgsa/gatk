## (howto) Apply hard filters to a call set

http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

<h4>Objective</h4>
<p>Apply hard filters to a variant callset that is too small for VQSR or for which truth/training sets are not available.</p>
<h4>Caveat</h4>
<p>This document is intended to illustrate how to compose and run the commands involved in applying the hard filtering method. The annotations and values used may not reflect the most recent recommendations. Be sure to read the documentation about <a href="https://www.broadinstitute.org/gatk/guide/article?id=3225">why you would use hard filters</a> and <a href="https://www.broadinstitute.org/gatk/guide/article?id=6925">how to understand and improve upon the generic hard filtering recommendations</a> that we provide.</p>
<h4>Steps</h4>
<ol>
<li>Extract the SNPs from the call set</li>
<li>Determine parameters for filtering SNPs</li>
<li>Apply the filter to the SNP call set </li>
<li>Extract the Indels from the call set</li>
<li>Determine parameters for filtering indels</li>
<li>Apply the filter to the Indel call set </li>
</ol>
<hr />
<h3>1. Extract the SNPs from the call set</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T SelectVariants \ 
    -R reference.fa \ 
    -V raw_variants.vcf \ 
    -selectType SNP \ 
    -o raw_snps.vcf </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a VCF file called <code>raw_snps.vcf</code>, containing just the SNPs from the original file of raw variants. </p>
<hr />
<h3>2. Determine parameters for filtering SNPs</h3>
<p>SNPs matching any of these conditions will be considered bad and filtered out, <em>i.e.</em> marked <code>FILTER</code> in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the SNP using the culprit annotation. SNPs that do not match any of these conditions will be considered good and marked <code>PASS</code> in the output VCF file. </p>
<ul>
<li>QualByDepth (QD) 2.0</li>
</ul>
<p>This is the variant confidence (from the <code>QUAL</code> field) divided by the unfiltered depth of non-reference samples.</p>
<ul>
<li>FisherStrand (FS) 60.0</li>
</ul>
<p>Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.</p>
<ul>
<li>RMSMappingQuality (MQ) 40.0</li>
</ul>
<p>This is the Root Mean Square of the mapping quality of the reads across all samples.</p>
<ul>
<li>MappingQualityRankSumTest (MQRankSum) -12.5</li>
</ul>
<p>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, <em>i.e.</em> this will only be applied to heterozygous calls. </p>
<ul>
<li>ReadPosRankSumTest (ReadPosRankSum) -8.0</li>
</ul>
<p>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, <em>i.e.</em> this will only be applied to heterozygous calls.</p>
<hr />
<h3>3. Apply the filter to the SNP call set</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T VariantFiltration \ 
    -R reference.fa \ 
    -V raw_snps.vcf \ 
    --filterExpression "QD &lt; 2.0 || FS &gt; 60.0 || MQ &lt; 40.0 || MQRankSum &lt; -12.5 || ReadPosRankSum &lt; -8.0" \ 
    --filterName "my_snp_filter" \ 
    -o filtered_snps.vcf </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a VCF file called <code>filtered_snps.vcf</code>, containing all the original SNPs from the <code>raw_snps.vcf</code> file, but now the SNPs are annotated with either <code>PASS</code> or <code>FILTER</code> depending on whether or not they passed the filters. </p>
<p>For SNPs that failed the filter, the variant annotation also includes the name of the filter. That way, if you apply several different filters (simultaneously or sequentially), you can keep track of which filter(s) each SNP failed, and later you can retrieve specific subsets of your calls using the SelectVariants tool. To learn more about composing different types of filtering expressions and retrieving subsets of variants using SelectVariants, please see the online GATK documentation. </p>
<hr />
<h3>4. Extract the Indels from the call set</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T SelectVariants \ 
    -R reference.fa \ 
    -V raw_HC_variants.vcf \ 
    -selectType INDEL \ 
    -o raw_indels.vcf </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a VCF file called <code>raw_indels.vcf</code>, containing just the Indels from the original file of raw variants.</p>
<hr />
<h3>5. Determine parameters for filtering Indels.</h3>
<p>Indels matching any of these conditions will be considered bad and filtered out, <em>i.e.</em> marked <code>FILTER</code> in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the indel using the culprit annotation. Indels that do not match any of these conditions will be considered good and marked <code>PASS</code> in the output VCF file. </p>
<ul>
<li>QualByDepth (QD) 2.0</li>
</ul>
<p>This is the variant confidence (from the <code>QUAL</code> field) divided by the unfiltered depth of non-reference samples.</p>
<ul>
<li>FisherStrand (FS) 200.0</li>
</ul>
<p>Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.</p>
<ul>
<li>ReadPosRankSumTest (ReadPosRankSum) 20.0</li>
</ul>
<p>This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, <em>i.e.</em> this will only be applied to heterozygous calls.</p>
<hr />
<h3>6. Apply the filter to the Indel call set</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T VariantFiltration \ 
    -R reference.fa \ 
    -V raw_indels.vcf \ 
    --filterExpression "QD &lt; 2.0 || FS &gt; 200.0 || ReadPosRankSum &lt; -20.0" \ 
    --filterName "my_indel_filter" \ 
    -o filtered_indels.vcf </code class="pre_md"></pre>
<h4>Expected Result</h4>
<p>This creates a VCF file called <code>filtered_indels.vcf</code>, containing all the original Indels from the <code>raw_indels.vcf</code> file, but now the Indels are annotated with either <code>PASS</code> or <code>FILTER</code> depending on whether or not they passed the filters.</p>
<p>For Indels that failed the filter, the variant annotation also includes the name of the filter. That way, if you apply several different filters (simultaneously or sequentially), you can keep track of which filter(s) each Indel failed, and later you can retrieve specific subsets of your calls using the SelectVariants tool. To learn more about composing different types of filtering expressions and retrieving subsets of variants using SelectVariants, please see the online GATK documentation. </p>