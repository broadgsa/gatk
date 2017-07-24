## What should I use as known variants/sites for running tool X?

http://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x

<h3>1. Notes on known sites</h3>
<h4>Why are they important?</h4>
<p>Each tool uses known sites differently, but what is common to all is that they use them to help distinguish true variants from false positives, which is very important to how these tools work. If you don't provide known sites, the statistical analysis of the data will be skewed, which can dramatically affect the sensitivity and reliability of the results. </p>
<p>In the variant calling pipeline, the only tools that do not strictly require known sites are UnifiedGenotyper and HaplotypeCaller.</p>
<h4>Human genomes</h4>
<p>If you're working on human genomes, you're in luck. We provide sets of known sites in the human genome as part of our <a href="http://www.broadinstitute.org/gatk/guide/article?id=1213">resource bundle</a>, and we can give you specific Best Practices recommendations on which sets to use for each tool in the variant calling pipeline. See the next section for details.</p>
<h4>Non-human genomes</h4>
<p>If you're working on genomes of other organisms, things may be a little harder -- but don't panic, we'll try to help as much as we can. We've started a community discussion in the forum on <a href="http://gatkforums.broadinstitute.org/discussion/1243">What are the standard resources for non-human genomes?</a> in which we hope people with non-human genomics experience will share their knowledge. </p>
<p>And if it turns out that there is as yet no suitable set of known sites for your organisms, here's how to make your own for the purposes of BaseRecalibration: First, do an initial round of SNP calling on your original, unrecalibrated data. Then take the SNPs that you have the highest confidence in and use that set as the database of known SNPs by feeding it as a VCF file to the base quality score recalibrator. Finally, do a real round of SNP calling with the recalibrated data. These steps could be repeated several times until convergence. Good luck!</p>
<p>Some experimentation will be required to figure out the best way to find the highest confidence SNPs for use here. Perhaps one could call variants with several different calling algorithms and take the set intersection. Or perhaps one could do a very strict round of filtering and take only those variants which pass the test. </p>
<h3>2. Recommended sets of known sites per tool</h3>
<h4>Summary table</h4>
<table class="table table-striped">
<thead>
<tr>
<th style="text-align: left;"><strong>Tool</strong></th>
<th style="text-align: center;"><strong>dbSNP 129</strong></th>
<th style="text-align: center;"><strong>dbSNP &gt;132</strong></th>
<th style="text-align: center;"><strong>Mills indels</strong></th>
<th style="text-align: center;"><strong>1KG indels</strong></th>
<th style="text-align: center;"><strong>HapMap</strong></th>
<th style="text-align: center;"><strong>Omni</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">RealignerTargetCreator</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr>
<td style="text-align: left;">IndelRealigner</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr>
<td style="text-align: left;">BaseRecalibrator</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr>
<td style="text-align: left;">(UnifiedGenotyper/ HaplotypeCaller)</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr>
<td style="text-align: left;">VariantRecalibrator</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;">X</td>
<td style="text-align: center;">X</td>
</tr>
<tr>
<td style="text-align: left;">VariantEval</td>
<td style="text-align: center;">X</td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
</tbody>
</table>
<h4>RealignerTargetCreator and IndelRealigner</h4>
<p>These tools require known indels passed with the <code>-known</code> argument to function properly. We use both the following files:</p>
<ul>
<li>Mills_and_1000G_gold_standard.indels.b37.sites.vcf</li>
<li>1000G_phase1.indels.b37.vcf (currently from the 1000 Genomes Phase I indel calls)</li>
</ul>
<h4>BaseRecalibrator</h4>
<p>This tool requires known SNPs and indels passed with the <code>-knownSites</code> argument to function properly. We use all the following files:</p>
<ul>
<li>The most recent dbSNP release (build ID &gt; 132)</li>
<li>Mills_and_1000G_gold_standard.indels.b37.sites.vcf</li>
<li>1000G_phase1.indels.b37.vcf (currently from the 1000 Genomes Phase I indel calls)</li>
</ul>
<h4>UnifiedGenotyper / HaplotypeCaller</h4>
<p>These tools do NOT require known sites, but if SNPs are provided with the <code>-dbsnp</code> argument they will use them for variant annotation. We use this file:</p>
<ul>
<li>The most recent dbSNP release (build ID &gt; 132)</li>
</ul>
<h4>VariantRecalibrator</h4>
<p>For VariantRecalibrator, please see the <a href="http://www.broadinstitute.org/gatk/guide/article?id=1259">FAQ article on VQSR training sets and arguments</a>.</p>
<h4>VariantEval</h4>
<p>This tool requires known SNPs passed with the <code>-dbsnp</code> argument to function properly. We use the following file:</p>
<ul>
<li>A version of dbSNP subsetted to only sites discovered in or before dbSNP BuildID 129, which excludes the impact of the 1000 Genomes project and is useful for evaluation of dbSNP rate and Ti/Tv values at novel sites.</li>
</ul>