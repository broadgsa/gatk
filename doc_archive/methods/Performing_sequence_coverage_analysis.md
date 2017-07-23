## Performing sequence coverage analysis

http://gatkforums.broadinstitute.org/gatk/discussion/40/performing-sequence-coverage-analysis

<h3>Overview</h3>
<p>This document describes the tools and concepts involved in performing sequence coverage analysis, where the purpose is to answer the common question: &quot;(Where) Do I have enough sequence data to be empowered to discover variants with reasonable confidence?&quot;. </p>
<p>The tools involved are the following:</p>
<ul>
<li>
<p><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php">DepthOfCoverage</a>:</strong> for QC'ing coverage in whole-genome data (WGS)</p>
</li>
<li><strong><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_diagnostics_diagnosetargets_DiagnoseTargets.php">DiagnoseTargets</a>:</strong> for QC'ing coverage in exome data (WEx) </li>
</ul>
<p>For an overview of the major annotations that are used by variant callers to express read depth at a variant site, and guidelines for using those metrics to evaluate variants, please see <a href="https://www.broadinstitute.org/gatk/guide/article?id=4721">this document</a>.</p>
<hr />
<h3>Introduction to coverage analysis as a QC method</h3>
<p>Coverage analysis generally aims to answer the common question: &quot;(Where) Do I have enough sequence data to be empowered to discover variants with reasonable confidence?&quot;. </p>
<p><strong>This section is incomplete.</strong></p>
<hr />
<h3>Using DepthOfCoverage to QC whole-genome data</h3>
<p><a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.html">DepthOfCoverage</a> is a coverage profiler for a (possibly multi-sample) bam file. It uses a granular histogram that can be user-specified to present useful aggregate coverage data. It reports the following metrics over the entire .bam file:</p>
<ul>
<li>Total, mean, median, and quartiles for each partition type: aggregate</li>
<li>Total, mean, median, and quartiles for each partition type: for each interval</li>
<li>A series of histograms of the number of bases covered to Y depth for each partition type (granular; e.g. Y can be a range, like 16 to 22)</li>
<li>A matrix of counts of the number of intervals for which at least Y samples and/or read groups had a median coverage of at least X</li>
<li>A matrix of counts of the number of bases that were covered to at least X depth, in at least Y groups (e.g. # of loci with ≥15x coverage for ≥12 samples)</li>
<li>A matrix of proportions of the number of bases that were covered to at least X depth, in at least Y groups (e.g. proportion of loci with ≥18x coverage for ≥15 libraries)</li>
</ul>
<p>That last matrix is key to answering the question posed above, so we recommend running this tool on all samples together.</p>
<p>Note that DepthOfCoverage can be configured to output these statistics aggregated over genes by providing it with a RefSeq gene list.</p>
<p>DepthOfCoverage also outputs, by default, the total coverage at every locus, and the coverage per sample and/or read group. This behavior can optionally be turned off, or switched to base count mode, where base counts will be output at each locus, rather than total depth.</p>
<p>To get a summary of coverage by each gene, you may supply a refseq (or alternative) gene list via the argument</p>
<pre><code class="pre_md">-geneList /path/to/gene/list.txt</code class="pre_md"></pre>
<p>The provided gene list must be of the following format:</p>
<pre><code class="pre_md">585     NM_001005484    chr1    +       58953   59871   58953   59871   1       58953,  59871,  0       OR4F5   cmpl    cmpl    0,
587     NM_001005224    chr1    +       357521  358460  357521  358460  1       357521, 358460, 0       OR4F3   cmpl    cmpl    0,
587     NM_001005277    chr1    +       357521  358460  357521  358460  1       357521, 358460, 0       OR4F16  cmpl    cmpl    0,
587     NM_001005221    chr1    +       357521  358460  357521  358460  1       357521, 358460, 0       OR4F29  cmpl    cmpl    0,
589     NM_001005224    chr1    -       610958  611897  610958  611897  1       610958, 611897, 0       OR4F3   cmpl    cmpl    0,
589     NM_001005277    chr1    -       610958  611897  610958  611897  1       610958, 611897, 0       OR4F16  cmpl    cmpl    0,
589     NM_001005221    chr1    -       610958  611897  610958  611897  1       610958, 611897, 0       OR4F29  cmpl    cmpl    0,</code class="pre_md"></pre>
<p>For users who have access to internal Broad resources, the properly-formatted file containing refseq genes and transcripts is located at</p>
<pre><code class="pre_md">/humgen/gsa-hpprojects/GATK/data/refGene.sorted.txt</code class="pre_md"></pre>
<p>If you do not have access (if you don't know, you probably don't have it), you can generate your own as described <a href="https://www.broadinstitute.org/gatk/guide/article?id=1329">here</a>.</p>
<p>If you supply the <code>-geneList</code> argument, DepthOfCoverage will output an additional summary file that looks as follows:</p>
<pre><code class="pre_md">Gene_Name     Total_Cvg       Avg_Cvg       Sample_1_Total_Cvg    Sample_1_Avg_Cvg    Sample_1_Cvg_Q3       Sample_1_Cvg_Median      Sample_1_Cvg_Q1
SORT1    594710  238.27  594710  238.27  165     245     330
NOTCH2  3011542 357.84  3011542 357.84  222     399     &amp;gt;500
LMNA    563183  186.73  563183  186.73  116     187     262
NOS1AP  513031  203.50  513031  203.50  91      191     290</code class="pre_md"></pre>
<p>Note that the gene coverage will be aggregated only over samples (not read groups, libraries, or other types). The <code>-geneList</code> argument also requires specific intervals within genes to be given (say, the particular exons you are interested in, or the entire gene), and it functions by aggregating coverage from the interval level to the gene level, by referencing each interval to the gene in which it falls. Because by-gene aggregation looks for intervals that overlap genes, <code>-geneList</code> is ignored if <code>-omitIntervals</code> is thrown.</p>
<hr />
<h3>Using DiagnoseTargets to QC whole-exome data</h3>
<p><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_diagnostics_diagnosetargets_DiagnoseTargets.php">DiagnoseTargets</a> produces a pseudo-VCF file that provides a &quot;CallableStatus&quot; judgment for each position or range of positions in the input bam file. The possible judgments are as follows:</p>
<ul>
<li>
<p>PASS : The base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE.</p>
</li>
<li>
<p>COVERAGE_GAPS : Absolutely no coverage was observed at a locus, regardless of the filtering parameters.</p>
</li>
<li>
<p>LOW_COVERAGE : There were less than min. depth bases at the locus, after applying filters.</p>
</li>
<li>
<p>EXCESSIVE_COVERAGE: More than <code>-maxDepth</code> read at the locus, indicating some sort of mapping problem.</p>
</li>
<li>
<p>POOR_QUALITY : More than <code>--maxFractionOfReadsWithLowMAPQ</code> at the locus, indicating a poor mapping quality of the reads.</p>
</li>
<li>
<p>BAD_MATE : The reads are not properly mated, suggesting mapping errors.</p>
</li>
<li>NO_READS : There are no reads contained in the interval.</li>
</ul>