## (howto) Call variants with the UnifiedGenotyper

http://gatkforums.broadinstitute.org/gatk/discussion/2804/howto-call-variants-with-the-unifiedgenotyper

<h3>Note: the UnifiedGenotyper has been replaced by HaplotypeCaller, which is a much better tool. UG is still available but you should really consider using HC instead.</h3>
<h4>Objective</h4>
<p>Call variants on a haploid genome with the UnifiedGenotyper, producing a raw (unfiltered) VCF.</p>
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
<li>Ploidy (<code>-ploidy</code>) </li>
</ul>
<p>In its basic use, this is the ploidy (number of chromosomes) per sample. By default it is set to 2, to process diploid organisms' genomes, but it can be set to any other desired value, starting at 1 for haploid organisms, and up for polyploids. This argument can also be used to handle pooled data. For that purpose, you'll need to set <code>-ploidy</code> to <code>Number of samples in each pool * Sample Ploidy</code>. There is no fixed upper limit, but keep in mind that high-level ploidy will increase processing times since the calculations involved are more complex. For full details on how to process pooled data, see Eran et al. (in preparation).</p>
<ul>
<li>Genotype likelihood model (<code>-glm</code>) </li>
</ul>
<p>This is the model that the program will use to calculate the genotype likelihoods. By default, it is set to <code>SNP</code>, but it can also be set to <code>INDEL</code> or <code>BOTH</code>. If set to <code>BOTH</code>, both SNPs and Indels will be called in the same run and be output to the same variants file.</p>
<ul>
<li>Emission confidence threshold (<code>-stand_emit_conf</code>) </li>
</ul>
<p>This is the minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant.</p>
<ul>
<li>Calling confidence threshold (<code>-stand_call_conf</code>) </li>
</ul>
<p>This is the minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. If a site's associated genotype has a confidence score lower than the calling threshold, the program will emit the site as filtered and will annotate it as LowQual. This threshold separates high confidence calls from low confidence calls.</p>
<p><em>The terms called and filtered are tricky because they can mean different things depending on context. In ordinary language, people often say a site was called if it was emitted as variant. But in the GATK's technical language, saying a site was called means that that site passed the confidence threshold test. For filtered, it's even more confusing, because in ordinary language, when people say that sites were filtered, they usually mean that those sites successfully passed a filtering test. However, in the GATK's technical language, the same phrase (saying that sites were filtered) means that those sites failed the filtering test. In effect, it means that those would be filtered out if the filter was used to actually remove low-confidence calls from the callset, instead of just tagging them. In both cases, both usages are valid depending on the point of view of the person who is reporting the results. So it's always important to check what is the context when interpreting results that include these terms.</em></p>
<hr />
<h3>2. Call variants in your sequence data</h3>
<p>Run the following GATK command: </p>
<pre><code class="pre_md">java -jar GenomeAnalysisTK.jar \ 
    -T UnifiedGenotyper \ 
    -R haploid_reference.fa \ 
    -I haploid_reads.bam \ 
    -L 20 \ 
    -glm BOTH \ 
    --stand_call_conf 30 \ 
    --stand_emit_conf 10 \ 
    -o raw_ug_variants.vcf </code class="pre_md"></pre>
<p>This creates a VCF file called <code>raw_ug_variants.vcf</code>, containing all the sites that the UnifiedGenotyper evaluated to be potentially variant. </p>
<p><em>Note that <code>-L</code> specifies that we only want to run the command on a subset of the data (here, chromosome 20). This is useful for testing as well as other purposes. For example, when running on exome data, we use <code>-L</code> to specify a file containing the list of exome targets corresponding to the capture kit that was used to generate the exome libraries.</em></p>
<p>Although you now have a nice fresh set of variant calls, the variant discovery stage is not over. The distinctions made by the caller itself between low-confidence calls and the rest is very primitive, and should not be taken as a definitive guide for filtering. The GATK callers are designed to be very lenient in calling variants, so it is extremely important to apply one of the recommended filtering methods (variant recalibration or hard-filtering), in order to move on to downstream analyses with the highest-quality call set possible.</p>