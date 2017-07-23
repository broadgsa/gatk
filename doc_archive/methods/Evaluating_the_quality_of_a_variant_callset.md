## Evaluating the quality of a variant callset

http://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset

<h2>Introduction</h2>
<p>Running through the steps involved in <a href="https://www.broadinstitute.org/gatk/guide/bp_step.php?p=2">variant discovery</a> (calling variants, joint genotyping and applying filters) produces a variant callset in the form of a VCF file. So what’s next? Technically, that callset is ready to be used in downstream analysis. But before you do that, we recommend running some quality control analyses to evaluate how “good” that callset is. </p>
<p>To be frank, distinguishing between a “good” callset and a “bad” callset is a complex problem. If you knew the absolute truth of what variants are present or not in your samples, you probably wouldn’t be here running variant discovery on some high-throughput sequencing data. Your fresh new callset is your attempt to discover that truth. So how do you know how close you got?</p>
<h3>Methods for variant evaluation</h3>
<p>There are several methods that you can apply which offer different insights into the probable biological truth, all with their own pros and cons. Possibly the most trusted method is Sanger sequencing of regions surrounding putative variants. However, it is also the least scalable as it would be prohibitively costly and time-consuming to apply to an entire callset. Typically, Sanger sequencing is only applied to validate candidate variants that are judged highly likely. Another popular method is to evaluate concordance against results obtained from a genotyping chip run on the same samples. This is much more scalable, and conveniently also doubles as a quality control method to detect sample swaps. Although it only covers the subset of known variants that the chip was designed for, this method can give you a pretty good indication of both sensitivity (ability to detect true variants) and specificity (not calling variants where there are none). This is something we do systematically for all samples in the Broad’s production pipelines.</p>
<p>The third method, presented here, is to evaluate how your variant callset stacks up against another variant callset (typically derived from other samples) that is considered to be a <strong>truth set</strong> (sometimes referred to as a <strong>gold standard</strong> -- these terms are very close and often used interchangeably). The general idea is that key properties of your callset (metrics discussed later in the text) should roughly match those of the truth set. This method is not meant to render any judgments about the veracity of individual variant calls; instead, it aims to estimate the overall quality of your callset and detect any red flags that might be indicative of error.</p>
<h3>Underlying assumptions and truthiness<sup>*</sup>: a note of caution</h3>
<p>It should be immediately obvious that there are two important assumptions being made here: <strong>1</strong>) that the content of the truth set has been validated somehow and is considered especially trustworthy; and <strong>2</strong>) that your samples are expected to have similar genomic content as the population of samples that was used to produce the truth set. These assumptions are not always well-supported, depending on the truth set, your callset, and what they have (or don’t have) in common. You should always keep this in mind when choosing a truth set for your evaluation; it’s a jungle out there. Consider that if anyone can submit variants to a truth set’s database without a well-regulated validation process, and there is no process for removing variants if someone later finds they were wrong (I’m looking at you, dbSNP), you should be extra cautious in interpreting results.
<sup>*With apologies to <a href="https://en.wikipedia.org/wiki/Truthiness">Stephen Colbert</a>.</sup></p>
<h3>Validation</h3>
<p>So what constitutes validation? Well, the best validation is done with orthogonal methods, meaning that it is done with technology (wetware, hardware, software, etc.) that is not subject to the same error modes as the sequencing process. Calling variants with two callers that use similar algorithms? Great way to reinforce your biases. It won’t mean anything that both give the same results; they could both be making the same mistakes. On the wetlab side, Sanger and genotyping chips are great validation tools; the technology is pretty different, so they tend to make different mistakes. Therefore it means more if they agree or disagree with calls made from high-throughput sequencing. </p>
<h3>Matching populations</h3>
<p>Regarding the population genomics aspect: it’s complicated -- especially if we’re talking about humans (I am). There’s a lot of interesting literature on this topic; for now let’s just summarize by saying that some important variant calling metrics vary depending on ethnicity. So if you are studying a population with a very specific ethnic composition, you should try to find a truth set composed of individuals with a similar ethnic background, and adjust your expectations accordingly for some metrics.</p>
<p>Similar principles apply to non-human genomic data, with important variations depending on whether you’re looking at wild or domesticated populations, natural or experimentally manipulated lineages, and so on. Unfortunately we can’t currently provide any detailed guidance on this topic, but hopefully this explanation of the logic and considerations involved will help you formulate a variant evaluation strategy that is appropriate for your organism of interest.</p>
<hr />
<h2>Variant evaluation metrics</h2>
<p>So let’s say you’ve got your fresh new callset and you’ve found an appropriate truth set. You’re ready to look at some metrics (but don’t worry yet about how; we’ll get to that soon enough). There are several metrics that we recommend examining in order to evaluate your data. The set described here should be considered a minimum and is by no means exclusive. It is nearly always better to evaluate more metrics if you possess the appropriate data to do so -- and as long as you understand why those additional metrics are meaningful. Please don’t try to use metrics that you don’t understand properly, because misunderstandings lead to confusion; confusion leads to worry; and worry leads to too many desperate posts on the GATK forum. </p>
<h3>Variant-level concordance and genotype concordance</h3>
<p>The relationship between variant-level concordance and genotype concordance is illustrated in <a href="https://us.v-cdn.net/5019796/uploads/FileUpload/09/6ba291fb1b8fe47895208d5e1bf380.png">this figure</a>.</p>
<ul>
<li>
<p><strong>Variant-level concordance</strong> (aka % Concordance) gives the percentage of variants in your samples that match (are concordant with) variants in your truth set. It essentially serves as a check of how well your analysis pipeline identified variants contained in the truth set. Depending on what you are evaluating and comparing, the interpretation of percent concordance can vary quite significantly.
Comparing your sample(s) against genotyping chip results matched per sample allows you to evaluate whether you missed any real variants within the scope of what is represented on the chip. Based on that concordance result, you can extrapolate what proportion you may have missed out of the real variants not represented on the chip.
If you don't have a sample-matched truth set and you're comparing your sample against a truth set derived from a population, your interpretation of percent concordance will be more limited. You have to account for the fact that some variants that are real in your sample will not be present in the population and that conversely, many variants that are in the population will not be present in your sample. In both cases, &quot;how many&quot; depends on how big the population is and how representative it is of your sample's background.
Keep in mind that for most tools that calculate this metric, all unmatched variants (present in your sample but not in the truth set) are considered to be false positives. Depending on your trust in the truth set and whether or not you expect to see true, novel variants, these unmatched variants could warrant further investigation -- or they could be artifacts that you should ignore.</p>
</li>
<li><strong>Genotype concordance</strong> is a similar metric but operates at the genotype level. It allows you to evaluate, within a set of variant calls that are present in both your sample callset and your truth set, what proportion of the genotype calls have been assigned correctly. This assumes that you are comparing your sample to a matched truth set derived from the same original sample. </li>
</ul>
<h3>Number of Indels &amp; SNPs and TiTv Ratio</h3>
<p>These metrics are widely applicable. The table below summarizes their expected value ranges for Human Germline Data:</p>
<table class="table table-striped">
<thead>
<tr>
<th><font size = "3"> Sequencing Type</font></th>
<th><font size = "3"># of Variants*</font></th>
<th><font size = "3">TiTv Ratio</font></th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>WGS</strong></td>
<td>~4.4M</td>
<td>2.0-2.1</td>
</tr>
<tr>
<td><strong>WES</strong></td>
<td>~41k</td>
<td>3.0-3.3</td>
</tr>
</tbody>
</table>
<p><sup>*for a single sample</sup></p>
<ul>
<li>
<p><strong>Number of Indels &amp; SNPs</strong>
The number of variants detected in your sample(s) are counted separately as indels (<strong>in</strong>sertions and <strong>del</strong>etions) and SNPs (<strong>S</strong>ingle <strong>N</strong>ucleotide <strong>P</strong>olymorphism<strong>s</strong>). Many factors can affect this statistic including whole exome (WES) versus whole genome (WGS) data, cohort size, strictness of filtering through the GATK pipeline, the ethnicity of your sample(s), and even algorithm improvement due to a software update. For reference, Nature's recently published <a href="http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html">2015 paper</a> in which various ethnicities in a moderately large cohort were analyzed for number of variants. As such, this metric alone is insufficient to confirm data validity, but it can raise warning flags when something went extremely wrong: e.g. 1000 variants in a large cohort WGS data set, or 4 billion variants in a ten-sample whole-exome set.</p>
</li>
<li><strong>TiTv Ratio</strong>
This metric is the ratio of <strong>t</strong>rans<strong>i</strong>tion (Ti) to <strong>t</strong>rans<strong>v</strong>ersion (Tv) SNPs. If the distribution of transition and transversion mutations were random (i.e. without any biological influence) we would expect a ratio of 0.5. This is simply due to the fact that there are twice as many possible transversion mutations than there are transitions. However, in the biological context, it is very common to see a methylated cytosine undergo deamination to become thymine. As this is a transition mutation, it has been shown to increase the expected random ratio from 0.5 to ~2.0<sup><a href="https://www.biostars.org/p/4751/">1</a></sup>. Furthermore, CpG islands, usually found in primer regions, have higher concentrations of methylcytosines. By including these regions, whole exome sequencing shows an even stronger lean towards transition mutations, with an expected ratio of 3.0-3.3. A significant deviation from the expected values could indicate artifactual variants causing bias. If your TiTv Ratio is too low, your callset likely has more false positives. <br>
It should also be noted that the TiTv ratio from exome-sequenced data will vary from the expected value based upon the length of flanking sequences. When we analyze exome sequence data, we add some padding (usually 100 bases) around the targeted regions (using the <code>-ip</code> engine argument) because this improves calling of variants that are at the edges of exons (whether inside the exon sequence or in the promoter/regulatory sequence before the exon). These flanking sequences are not subject to the same evolutionary pressures as the exons themselves, so the number of transition and transversion mutants lean away from the expected ratio. The amount of &quot;lean&quot; depends on how long the flanking sequence is.</li>
</ul>
<h3>Ratio of Insertions to Deletions (Indel Ratio)</h3>
<p>This metric is generally evaluated after filtering for purposes that are specific to your study, and the expected value range depends on whether you're looking for rare or common variants, as summarized in the table below.</p>
<table class="table table-striped">
<thead>
<tr>
<th><font size="3">Filtering for</font></th>
<th><font size="3">Indel Ratio</font></th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>common</strong></td>
<td>~1</td>
</tr>
<tr>
<td><strong>rare</strong></td>
<td>0.2-0.5</td>
</tr>
</tbody>
</table>
<p>A significant deviation from the expected ratios listed in the table above could indicate a bias resulting from artifactual variants.</p>
<hr />
<h2>Tools for performing variant evaluation</h2>
<h3><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php">VariantEval</a></h3>
<p>This is the GATK’s main tool for variant evaluation. It is designed to collect and calculate a variety of callset metrics that are organized in <strong>evaluation modules</strong>, which are listed in the <a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php">tool doc</a>. For each evaluation module that is enabled, the tool will produce a table containing the corresponding callset metrics based on the specified inputs (your callset of interest and one or more truth sets). By default, VariantEval will run with a specific subset of the available modules (listed below), but all evaluation modules can be enabled or disabled from the command line. We recommend setting the tool to produce only the metrics that you are interested in, because each active module adds to the computational requirements and overall runtime of the tool. </p>
<p>It should be noted that all module calculations only include variants that passed filtering (i.e. FILTER column in your vcf file should read PASS); variants tagged as filtered out will be ignored. It is not possible to modify this behavior. See the <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6211/howto-evaluate-a-callset-with-varianteval#latest">example analysis</a> for more details on how to use this tool and interpret its output.</p>
<h3><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php">GenotypeConcordance</a></h3>
<p>This tool calculates -- you’ve guessed it -- the genotype concordance between callsets. In earlier versions of GATK, GenotypeConcordance was itself a module within VariantEval. It was converted into a standalone tool to enable more complex genotype concordance calculations.</p>
<h3><a href="https://broadinstitute.github.io/picard/index.html">Picard tools</a></h3>
<p>The Picard toolkit includes two tools that perform similar functions to VariantEval and GenotypeConcordance, respectively called <a href="https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectVariantCallingMetrics.VariantCallingSummaryMetrics">CollectVariantCallingMetrics</a> and <a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html#GenotypeConcordanceDetailMetrics">GenotypeConcordance</a>. Both are relatively lightweight in comparison to their GATK equivalents; their functionalities are more limited, but they do run quite a bit faster. See the <a href="http://gatkforums.broadinstitute.org/gatk/discussion/6186/howto-evaluate-a-callset-with-collectvariantcallingmetrics#latest">example analysis</a> of CollectVariantCallingMetrics for details on its use and data interpretation. Note that in the coming months, the Picard tools are going to be integrated into the next major version of GATK, so at that occasion we plan to consolidate these two pairs of homologous tools to eliminate redundancy.</p>
<h3>Which tool should I use?</h3>
<p>We recommend Picard's version of each tool for most cases. The GenotypeConcordance tools provide mostly the same information, but Picard's version is preferred by Broadies. Both VariantEval and CollectVariantCallingMetrics produce similar metrics, however the latter runs faster and is scales better for larger cohorts. By default, CollectVariantCallingMetrics stratifies by sample, allowing you to see the value of relevant statistics as they pertain to specific samples in your cohort. It includes all metrics discussed here, as well as a few more. On the other hand, VariantEval provides many more metrics beyond the minimum described here for analysis. It should be noted that none of these tools use phasing to determine metrics. </p>
<p><strong>So when should I use CollectVariantCallingMetrics?</strong></p>
<ul>
<li>If you have a very large callset</li>
<li>If you want to look at the metrics discussed here and not much else</li>
<li>If you want your analysis back quickly</li>
</ul>
<p><strong>When should I use VariantEval?</strong></p>
<ul>
<li>When you require a more detailed analysis of your callset</li>
<li>If you need to stratify your callset by another factor (allele frequency, indel size, etc.)</li>
<li>If you need to compare to multiple truth sets at the same time</li>
</ul>