## How to use Firepony

http://gatkforums.broadinstitute.org/gatk/discussion/6021/how-to-use-firepony

<p>Firepony can be run with the following command line arguments:</p>
<pre><code class="pre_md">firepony -r &lt;reference FASTA file&gt; -s &lt;SNP database file&gt; -o &lt;output table file&gt; &lt;input alignment file&gt;</code class="pre_md"></pre>
<p>where:</p>
<ul>
<li><code>-r</code> specifies the path to the reference file (in uncompressed FASTA format, equivalent to GATK option <code>-R</code>)</li>
<li><code>-s</code> specifies the path to the SNP database file (in BCF or VCF format, equivalent to GATK option <code>-knownSites</code>). </li>
</ul>
<p>Firepony will load an index for the reference file if it exists, which enables on-demand loading of reference sequences as the SNP database is loaded.</p>
<p>For example, the following GATK command line:</p>
<pre><code class="pre_md">java -Xmx8g GenomeAnalysisTK-3.4.jar \
    -T BaseRecalibrator \
    -I NA12878D_HiSeqX_R1.deduplicated.bam \
    -R /store/ref/hs37d5.fa \
    -knownSites /store/dbsnp/dbsnp_138.b37.vcf \
    -o recal_data.table</code class="pre_md"></pre>
<p>would be replaced by the following Firepony command line:</p>
<pre><code class="pre_md">firepony \
    -r /store/ref/hs37d5.fa -s /store/dbsnp/dbsnp_138.b37.vcf \
    -o recal_data.table NA12878D_HiSeqX_R1.deduplicated.bam</code class="pre_md"></pre>
<p>Additional command line options are described in the help output for firepony invoked by </p>
<pre><code class="pre_md">`firepony --help`</code class="pre_md"></pre>
<p>Note that it is recommended to use the BCF format rather than VCF for SNP databases when running Firepony. Both generate the same results, but loading BCF files is much more efficient.</p>
<p>At the moment, Firepony only supports recalibrating Illumina reads with the default GATK BQSR parameters, listed below in BQSR table format. Expanding the parameter set as well as the number of supported instruments will be done based on user feedback.</p>
<pre><code class="pre_md">#:GATKTable:Arguments:Recalibration argument collection values used in this run
Argument                    Value
binary_tag_name             null
covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate
default_platform            null
deletions_default_quality   45
force_platform              null
indels_context_size         3
insertions_default_quality  45
low_quality_tail            2
maximum_cycle_value         500
mismatches_context_size     2
mismatches_default_quality  -1
no_standard_covs            false
quantizing_levels           16
recalibration_report        null
run_without_dbsnp           false
solid_nocall_strategy       THROW_EXCEPTION
solid_recal_mode            SET_Q_ZERO</code class="pre_md"></pre>