## Can I use GATK on non-diploid organisms?

http://gatkforums.broadinstitute.org/gatk/discussion/1214/can-i-use-gatk-on-non-diploid-organisms

<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/d3/424549fd16f54f89339a95b6634461.jpg" />
<p>In general most GATK tools don't care about ploidy. The major exception is, of course, at the variant calling step: the variant callers need to know what ploidy is assumed for a given sample in order to perform the appropriate calculations. </p>
<h3>Ploidy-related functionalities</h3>
<p>As of version 3.3, the HaplotypeCaller and GenotypeGVCFs are able to deal with non-diploid organisms (whether haploid or exotically polyploid). In the case of HaplotypeCaller, you need to specify the ploidy of your non-diploid sample with the <code>-ploidy</code> argument. HC can only deal with one ploidy at a time, so if you want to process different chromosomes with different ploidies (e.g. to call X and Y in males) you need to run them separately. On the bright side, you can combine the resulting files afterward. In particular, if you’re running the -ERC GVCF workflow, you’ll find that both CombineGVCFs and GenotypeGVCFs are able to handle mixed ploidies (between locations and between samples). Both tools are able to correctly work out the ploidy of any given sample at a given site based on the composition of the GT field, so they don’t require you to specify the -ploidy argument.</p>
<p>For earlier versions (all the way to 2.0) the fallback option is UnifiedGenotyper, which also accepts the <code>-ploidy</code> argument. </p>
<h3>Cases where ploidy needs to be specified</h3>
<ol>
<li>Native variant calling in haploid or polyploid organisms.  </li>
<li>Pooled calling where many pooled organisms share a single barcode and hence are treated as a single &quot;sample&quot;.  </li>
<li>Pooled validation/genotyping at known sites.  </li>
</ol>
<p>For normal organism ploidy, you just set the <code>-ploidy</code> argument to the desired number of chromosomes per organism. In the case of pooled sequencing experiments, this argument should be set to the number of chromosomes per barcoded sample, i.e. <code>(Ploidy per individual) * (Individuals in pool)</code>.</p>
<h2>Important limitations</h2>
<p>Several variant annotations are not appropriate for use with non-diploid cases. In particular, InbreedingCoeff will not be annotated on non-diploid calls. Annotations that do work and are supported in non-diploid use cases are the following: <code>QUAL</code>, <code>QD</code>, <code>SB</code>, <code>FS</code>, <code>AC</code>, <code>AF</code>, and Genotype annotations such as <code>PL</code>, <code>AD</code>, <code>GT</code>, etc.</p>
<p>You should also be aware of the fundamental accuracy limitations of high ploidy calling. Calling low-frequency variants in a pool or in an organism with high ploidy is hard because these rare variants become almost indistinguishable from sequencing errors. </p>