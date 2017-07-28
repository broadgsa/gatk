## I do not get the annotations I specified with -A

http://gatkforums.broadinstitute.org/gatk/discussion/6022/i-do-not-get-the-annotations-i-specified-with-a

<h3>The problem</h3>
<p>You specified <code>-A &lt;some annotation&gt;</code> in a command line invoking one of the annotation-capable tools (HaplotypeCaller, MuTect2, UnifiedGenotyper and VariantAnnotator), but that annotation did not show up in your output VCF. </p>
<p><em>Keep in mind that all annotations that are necessary to run our Best Practices are annotated by default, so you should generally not need to request annotations unless you're doing something a bit special.</em></p>
<h3>Why this happens &amp; solutions</h3>
<p>There can be several reasons why this happens, depending on the tool, the annotation, and you data. These are the four we see most often; if you encounter another that is not listed here, let us know in the comments.</p>
<ol>
<li>
<h4>You requested an annotation that cannot be calculated by the tool</h4>
<p>For example, you're running MuTect2 but requested an annotation that is specific to HaplotypeCaller. There should be an error message to that effect in the output log. It's not possible to override this; but if you believe the annotation should be available to the tool, let us know in the forum and we'll consider putting in a feature request. </p>
</li>
<li>
<h4>You requested an annotation that can only be calculated if an optional input is provided</h4>
<p>For example, you're running HaplotypeCaller and you want InbreedingCoefficient, but you didn't specify a pedigree file. There should be an error message to that effect in the output log. The solution is simply to provide the missing input file. Another example: you're running VariantAnnotator and you want to annotate Coverage, but you didn't specify a BAM file. The tool needs to see the read data in order to calculate the annotation, so again, you simply need to provide the BAM file. </p>
</li>
<li>
<h4>You requested an annotation that has requirements which are not met by some or all sites</h4>
<p>For example, you're looking at RankSumTest annotations, which require heterozygous sites in order to perform the necessary calculations, but you're running on haploid data so you don't have any het sites. There is no workaround; the annotation is not applicable to your data. Another example: you requested InbreedingCoefficient, but your population includes fewer than 10 founder samples, which are required for the annotation calculation. There is no workaround; the annotation is not applicable to your data.</p>
</li>
<li>
<h4>You requested an annotation that is already applied by default by the tool you are running</h4>
<p>For example, you requested Coverage from HaplotypeCaller, which already annotates this by default. There is currently a bug that causes some default annotations to be dropped from the list if specified on the command line. This will be addressed in an upcoming version. For now the workaround is to check what annotations are applied by default and NOT request them with <code>-A</code>. </p>
</li>
</ol>