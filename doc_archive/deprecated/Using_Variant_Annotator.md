## Using Variant Annotator

http://gatkforums.broadinstitute.org/gatk/discussion/49/using-variant-annotator

<h3>This document is out of date and has been retired. Please see the Annotation documentation in the Tool Docs as well as various other Guide articles for better materials on annotating variants.</h3>
<hr />
<p>2 SNPs with significant strand bias</p>
<img src="http://www.broadinstitute.org/gatk/media/pics/StrandFailure.png" />  
<p>Several SNPs with excessive coverage</p>
<img src="http://www.broadinstitute.org/gatk/media/pics/DoCFailure.png" />  
<p><strong>For a complete, detailed argument reference, refer to the GATK document page <a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_VariantAnnotator.html">here</a>.</strong></p>
<h3>Introduction</h3>
<p>In addition to true variation, variant callers emit a number of false-positives.  Some of these false-positives can be detected and rejected by various statistical tests.  VariantAnnotator provides a way of annotating variant calls as preparation for executing these tests.</p>
<p>Description of the haplotype score annotation</p>
<img src="http://www.broadinstitute.org/gatk/media/pics/HaplotypeScore.png" />  
<h3>Examples of Available Annotations</h3>
<p>The list below is not comprehensive.  Please use the <code>--list</code> argument to get a list of all possible annotations available.  Also, see <a href="http://www.broadinstitute.org/gatk/guide/article?id=1268">the FAQ article on understanding the Unified Genotyper's VCF files</a> for a description of some of the more standard annotations.</p>
<ul>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_BaseQualityRankSumTest.html">BaseQualityRankSumTest</a> (BaseQRankSum)</li>
<li><a href="https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">DepthOfCoverage</a> (DP)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_FisherStrand.html">FisherStrand</a> (FS)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_HaplotypeScore.html">HaplotypeScore</a> (HaplotypeScore)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityRankSumTest.html">MappingQualityRankSumTest</a> (MQRankSum)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_MappingQualityZero.html">MappingQualityZero</a> (MQ0)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_QualByDepth.html">QualByDepth</a> (QD)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_ReadPosRankSumTest.html">ReadPositionRankSumTest</a> (ReadPosRankSum)</li>
<li><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_RMSMappingQuality.html">RMSMappingQuality</a> (MQ)</li>
<li><a href="http://www.broadinstitute.org/gatk/guide/article?id=50">SnpEff</a>: Add genomic annotations using the third-party tool SnpEff with VariantAnnotator</li>
</ul>
<p>Note that technically the VariantAnnotator does not require reads (from a BAM file) to run; if no reads are provided, only those Annotations which don't use reads (e.g. Chromosome Counts) will be added. But most Annotations do require reads.  <strong>When running the tool we recommend that you add the <code>-L</code> argument with the variant rod to your command line for efficiency and speed.</strong></p>