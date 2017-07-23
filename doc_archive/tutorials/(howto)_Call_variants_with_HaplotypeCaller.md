## (howto) Call variants with HaplotypeCaller

http://gatkforums.broadinstitute.org/gatk/discussion/2803/howto-call-variants-with-haplotypecaller

<h4>Objective</h4>
<p>Call variants on a single genome with the HaplotypeCaller, producing a raw (unfiltered) VCF. </p>
<h4>Caveat</h4>
<p>This is meant only for single-sample analysis. To analyze multiple samples, see the Best Practices documentation on joint analysis.</p>
<h4>Prerequisites</h4>
<ul>
<li>TBD</li>
</ul>
<h4>Steps</h4>
<ol>
<li>Determine the basic parameters of the analysis</li>
<li>Call variants in your sequence data</li>
</ol>
<hr />
<h3>1. Determine the basic parameters of the analysis</h3>
<p>If you do not specify these parameters yourself, the program will use default values. However we recommend that you set them explicitly because it will help you understand how the results are bounded and how you can modify the program's behavior. </p>
<ul>
<li>Genotyping mode (<code>--genotyping_mode</code>) </li>
</ul>
<p>This specifies how we want the program to determine the alternate alleles to use for genotyping. In the default <code>DISCOVERY</code> mode, the program will choose the most likely alleles out of those it sees in the data. In <code>GENOTYPE_GIVEN_ALLELES</code> mode, the program will only use the alleles passed in from a VCF file (using the <code>-alleles</code> argument). This is useful if you just want to determine if a sample has a specific genotype of interest and you are not interested in other alleles. </p>
<ul>
<li>Emission confidence threshold (<code>-stand_emit_conf</code>) </li>
</ul>
<p>This is the minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant.</p>
<ul>
<li>Calling confidence threshold (<code>-stand_call_conf</code>) </li>
</ul>
<p>This is the minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. If a site's associated genotype has a confidence score lower than the calling threshold, the program will emit the site as filtered and will annotate it as LowQual. This threshold separates high confidence calls from low confidence calls.</p>
<p><em>The terms &quot;called&quot; and &quot;filtered&quot; are tricky because they can mean different things depending on context. In ordinary language, people often say a site was called if it was emitted as variant. But in the GATK's technical language, saying a site was called means that that site passed the confidence threshold test. For filtered, it's even more confusing, because in ordinary language, when people say that sites were filtered, they usually mean that those sites successfully passed a filtering test. However, in the GATK's technical language, the same phrase (saying that sites were filtered) means that those sites failed the filtering test. In effect, it means that those would be filtered out if the filter was used to actually remove low-confidence calls from the callset, instead of just tagging them. In both cases, both usages are valid depending on the point of view of the person who is reporting the results. So it's always important to check what is the context when interpreting results that include these terms.</em></p>
<hr />
<h3>2. Call variants in your sequence data</h3>
<h4>Action</h4>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T HaplotypeCaller \ 
    -R reference.fa \ 
    -I preprocessed_reads.bam \  
    -L 20 \ 
    --genotyping_mode DISCOVERY \ 
    -stand_emit_conf 10 \ 
    -stand_call_conf 30 \ 
    -o raw_variants.vcf </code class="pre_md"></pre>
<p><em>Note that <code>-L</code> specifies that we only want to run the command on a subset of the data (here, chromosome 20). This is useful for testing as well as other purposes, as documented <a href="http://www.broadinstitute.org/gatk/guide/article?id=4133">here</a>. For example, when running on exome data, we use <code>-L</code> to specify a file containing the list of exome targets corresponding to the capture kit that was used to generate the exome libraries.</em></p>
<h4>Expected Result</h4>
<p>This creates a VCF file called <code>raw_variants.vcf</code>, containing all the sites that the HaplotypeCaller evaluated to be potentially variant. Note that this file contains both SNPs and Indels.</p>
<p>Although you now have a nice fresh set of variant calls, the variant discovery stage is not over. The distinctions made by the caller itself between low-confidence calls and the rest is very primitive, and should not be taken as a definitive guide for filtering. The GATK callers are designed to be very lenient in calling variants, so it is extremely important to apply one of the recommended filtering methods (variant recalibration or hard-filtering), in order to move on to downstream analyses with the highest-quality call set possible.</p>