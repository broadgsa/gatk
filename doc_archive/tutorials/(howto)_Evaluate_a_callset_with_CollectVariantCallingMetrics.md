## (howto) Evaluate a callset with CollectVariantCallingMetrics

http://gatkforums.broadinstitute.org/gatk/discussion/6186/howto-evaluate-a-callset-with-collectvariantcallingmetrics

<h2>Related Documents</h2>
<ul>
<li><a href="https://www.broadinstitute.org/gatk/guide/article?id=6308">Evaluating the quality of a variant callset</a></li>
<li><a href="https://www.broadinstitute.org/gatk/guide/article?id=6211">(howto) Evaluate a callset with VariantEval</a></li>
</ul>
<h2>Context</h2>
<p>This document will walk you through use of Picard's CollectVariantCallingMetrics tool, an excellent tool for large callsets, especially if you need your results quickly and do not require many additional metrics to those described here. Your callset consists of variants identified by earlier steps in the <a href="https://www.broadinstitute.org/gatk/guide/best-practices">GATK best practices pipeline</a>, and now requires additional evaluation to determine where your callset falls on the spectrum of &quot;perfectly identifies all true, biological variants&quot; to &quot;only identifies artifactual or otherwise unreal variants&quot;. When variant calling, we want the callset to maximize the correct calls, while minimizing false positive calls. While very robust methods, such as Sanger sequencing, can be used to individually sequence each potential variant, statistical analysis can be used to evaluate callsets instead, saving both time and money. These callset-based analyses are accomplished by comparing relevant metrics between your samples and a known truth set, such as dbSNP. Two tools exist to examine these metrics: <a href="https://www.broadinstitute.org/gatk/guide/article?id=6211">VariantEval</a> in GATK, and CollectVariantCallingMetrics in Picard. While the latter is currently used in the Broad Institute's production pipeline, the merits to each tool, as well as the basis for variant evaluation, are discussed <a href="https://www.broadinstitute.org/gatk/guide/article?id=6308">here</a>. </p>
<hr />
<h2>Example Use</h2>
<h3>Command</h3>
<pre><code>java -jar picard.jar CollectVariantCallingMetrics \
INPUT=CEUtrio.vcf \
OUTPUT=CEUtrioMetrics \
DBSNP=dbsnp_138.b37.excluding_sites_after_129.vcf </code></pre>
<ul>
<li><strong>INPUT</strong>
The CEU trio (NA12892, NA12891, and 12878) from the 1000 Genomes Project is the input chosen for this example. It is the callset that we wish to examine the metrics on, and thus this is the field where you would specify the .vcf file containing your sample(s)'s variant calls.</li>
<li><strong>OUTPUT</strong>
The output for this command will be written to two files named CEUtrioMetrics.variant_calling_summary_metrics and CEUtrioMetrics.variant_calling_detail_metrics, hereafter referred to as &quot;summary&quot; and &quot;detail&quot;, respectively. The specification given in this field is applied as the name of the out files; the file extensions are provided by the tool itself.</li>
<li><strong>DBSNP</strong>
The last required input to run this tool is a dbSNP file. The one used here is available in the current <a href="https://www.broadinstitute.org/gatk/guide/article?id=1213">GATK bundle</a>. CollectVariantCallingMetrics utilizes this dbSNP file as a base of comparison against the sample(s) present in your vcf.</li>
</ul>
<h3>Getting Results</h3>
<p>After running the command, CollectVariantCallingMetrics will return both a detail and a summary metrics file. These files can be viewed as a text file if needed, or they can be read in as a table using your preferred spreadsheet viewer (e.g. Excel) or scripting language of your choice (e.g. python, R, etc.) The files contain headers and are tab-delimited; the commands for reading in the tables in RStudio are found below. (Note: Replace &quot;~/path/to/&quot; with the path to your output files as needed.)</p>
<pre><code>summary &lt;- read.table("~/path/to/CEUtrioMetrics.variant_calling_summary_metrics", header=TRUE, sep="\t")
detail &lt;- read.table("~/path/to/CEUtrioMetrics.variant_calling_detail_metrics", header=TRUE, sep="\t")</code></pre>
<ul>
<li><strong>Summary</strong>
The summary metrics file will contain a single row of data for each <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectVariantCallingMetrics.VariantCallingSummaryMetrics">metric</a>, taking into account all samples present in your <code>INPUT</code> file. </li>
<li><strong>Detail</strong>
The detail metrics file gives a breakdown of each <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectVariantCallingMetrics.VariantCallingSummaryMetrics">statistic</a> by sample. In addition to all metrics covered in the summary table, the detail table also contains entries for <code>SAMPLE_ALIAS</code> and <code>HET_HOMVAR_RATIO</code>. In the example case here, the detail file will contain metrics for the three different samples, NA12892, NA12891, and NA12878. </li>
</ul>
<h3>Analyzing Results</h3>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/82/07ef0007f402da4d56eba2ce50d9bb.png" />
<p><sup>*Concatenated in the above table are the detail file's (rows 1-3) and the summary file's (row 4) relevant metrics; for full output table, see attached image file.</sup></p>
<ul>
<li><strong>Number of Indels &amp; SNPs</strong>
This tool collects the number of SNPs (single nucleotide polymorphisms) and indels (insertions and deletions) as found in the variants file. It counts only biallelic sites and filters out multiallelic sites. Many factors affect these counts, including cohort size, relatedness between samples, strictness of filtering, ethnicity of samples, and even algorithm improvement due to updated software. While this metric alone is insufficient to evaluate your variants, it does provide a good baseline. It is reassuring to see that across the three related samples, we saw very similar numbers of SNPs and indels. It could be cause for concern if a particular sample had significantly more or fewer variants than the rest.</li>
<li><strong>Indel Ratio</strong>
The indel ratio is determined to be the total number of insertions divided by the total number of deletions; this tool does not include filtered variants in this calculation. Usually, the indel ratio is around 1, as insertions occur typically as frequently as deletions. However, in rare variant studies, indel ratio should be around 0.2-0.5. Our samples have an indel ratio of ~0.95, indicating that these variants are not likely to have a bias affecting their insertion/deletion ratio.</li>
<li><strong>TiTv Ratio</strong>
This metric is the ratio of transition (Ti) to transversion (Tv) mutations. For whole genome sequencing data, TiTv should be ~2.0-2.1, whereas whole exome sequencing data will have a TiTv ratio of ~3.0-3.3<sup><a href="http://www.nature.com/ng/journal/v43/n5/full/ng.806.html">1</a></sup>. In the case of the CEU trio of samples, the TiTv of ~2.08 and ~1.91 are within reason, and this variant callset is unlikely to have a bias affecting its transition/transversion ratio.</li>
</ul>