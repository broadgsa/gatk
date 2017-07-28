## Reference implementation: PairedEndSingleSampleWf pipeline

http://gatkforums.broadinstitute.org/gatk/discussion/7899/reference-implementation-pairedendsinglesamplewf-pipeline

<p><a name="top"></a></p>
<hr />
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/4b/1a77e7233c3f96dd41605eee2fa00a.png" align="right" height="270" style="margin:0px 0px 10px 5px"/> This document describes the workflow and details pertinent parameters of the <strong>PairedEndSingleSampleWf</strong> pipeline, which implements GATK Best Practices (ca. June 2016) for pre-processing human germline whole-genome sequencing (WGS) data. This pipeline uses GRCh38 as the reference genome and, as the name implies, is specific to processing paired end reads for a single sample. It begins with unaligned paired reads in BAM format and results in a sample-level SNP and INDEL variant callset in GVCF format. </p>
<ul>
<li>This document is specific to the public WDL script <strong>PublicPairedSingleSampleWf_160720.wdl</strong> with a July 20, 2016 date stamp found <a href="https://github.com/broadinstitute/wdl/blob/develop/scripts/broad_pipelines">here</a>.</li>
<li>The outline uses Docker container <strong>broadinstitute/genomes-in-the-cloud:2.2.2-1466113830</strong>. <a href="http://gatkforums.broadinstitute.org/wdl/discussion/6886">Docker</a> containers provide reproducible analysis environments, including specific tool versions. Thus, the commands and the given parameters pertain only to these specific tool versions. It is possible for prior or future versions of tools to give different results. </li>
<li>Furthermore, the parameters within the pipeline are optimized for WGS and GRCh38. Optimal performance on other types of data and other species may require changes to parameters. For example, <strong>Step 3</strong> includes code that calculates how to split the reference contigs for parallel processing that caps the genomic interval per process by the longest contig in the reference. This splits human GRCh38 18-ways and caps each group of contigs by the size of human chromosome 1. This approach may need revision for efficient processing of genomes where the longest contig length deviates greatly from human. </li>
</ul>
<p>The diagram above shows the relationship between the WORKFLOW steps that call on specific TASKS. Certain steps use genomic intervals to parallelize processes, and these are boxed in the workflow diagram. An overview of the data transformations is given in the WORKFLOW definitions section and granular details are given in the TASK definitions section in the order shown below. </p>
<hr />
<h2>Jump to a section</h2>
<h3>WORKFLOW definition <a href="#0">overview</a></h3>
<ol>
<li><a href="#1">Map with BWA-MEM and merge to create clean BAM</a></li>
<li><a href="#2">Flag duplicates with MarkDuplicates</a></li>
<li><a href="#3">Base quality score recalibration</a></li>
<li><a href="#4">Call SNP and INDEL variants with HaplotypeCaller</a></li>
</ol>
<h3>TASK definitions <a href="#tasksoverview">overview</a></h3>
<ul>
<li><a href="#GetBwaVersion">GetBwaVersion</a></li>
<li><a href="#SamToFastqAndBwaMem">SamToFastqAndBwaMem</a></li>
<li><a href="#MergeBamAlignment">MergeBamAlignment</a></li>
<li><a href="#MarkDuplicates">MarkDuplicates</a></li>
<li><a href="#SortAndFixTags">SortAndFixTags</a></li>
<li><a href="#CreateSequenceGroupingTSV">CreateSequenceGroupingTSV</a></li>
<li><a href="#BaseRecalibrator">BaseRecalibrator</a></li>
<li><a href="#GatherBqsrReports">GatherBqsrReports</a></li>
<li><a href="#ApplyBQSR">ApplyBQSR</a></li>
<li><a href="#GatherBamFiles">GatherBamFiles</a></li>
<li><a href="#ConvertToCram">ConvertToCram</a></li>
<li><a href="#HaplotypeCaller">HaplotypeCaller</a></li>
<li><a href="#GatherVCFs">GatherVCFs</a></li>
</ul>
<hr />
<h2>What is NOT covered</h2>
<ul>
<li>This document assumes a basic understanding of WDL components.</li>
<li>The JSON files describing inputs and outputs.</li>
<li>Runtime parameters optimized for Broad's Google Cloud Platform implementation. </li>
</ul>
<h3>Related resources</h3>
<ul>
<li>For details on interpreting and writing WDL scripts, see the <a href="https://software.broadinstitute.org/wdl/userguide/index">QuickStart guide</a>.</li>
<li>Scatter-Gather Parallelism. See <a href="https://software.broadinstitute.org/wdl/userguide/topic?name=wdl-plumbing">wdl plumbing options</a> for information.</li>
<li>Intervals lists. See <strong>Section 3</strong> of <a href="https://software.broadinstitute.org/gatk/documentation/article?id=1204">Article#1204</a>.</li>
</ul>
<hr />
<h2>Requirements</h2>
<h3>Software</h3>
<ul>
<li>See <a href="https://software.broadinstitute.org/wdl/guide/article?id=6671">Article#6671</a> for setup options in using Cromwell to execute WDL workflows. For an introduction, see <a href="https://www.broadinstitute.org/gatk/blog?id=7349">Article#7349</a>.</li>
<li>Docker container <code>broadinstitute/genomes-in-the-cloud:2.2.3-1469027018</code> uses the following tool versions for this pipeline. These tools in turn require Java JDK v8 (specifically 8u91) and Python v2.7.</li>
</ul>
<pre><code>DOCKER_VERSION="1.8.1"
PICARD_VERSION="1.1099"
GATK35_VERSION="3.5-0-g36282e4"
GATK4_VERSION="4.alpha-249-g7df4044"
SAMTOOLS_VER="1.3.1"
BWA_VER="0.7.13-r1126"</code></pre>
<h3>Scripts and data</h3>
<ul>
<li>Pipeline WDL script and JSON file defining input data.</li>
<li>(Optional) JSON file defining additional Cromwell options. </li>
<li>Human whole genome paired end sequence reads in unmapped BAM (uBAM) format. 
<ul>
<li>Each uBAM is per read group and the header defines read group <code>ID</code>, <code>SM</code>, <code>LB</code>, <code>PL</code> and optionally <code>PU</code>. Because each file is for the same sample, their <code>SM</code> fields will be identical. Each read has an <code>RG</code> tag.</li>
<li>Each uBAM has the same suffix, e.g. <code>.unmapped.bam</code>.</li>
<li>Each uBAM passes validation by <a href="https://software.broadinstitute.org/gatk/guide/article?id=7571">ValidateSamFile</a>.</li>
<li>Reads are in query-sorted order.</li>
</ul></li>
<li>GRCh38 reference genome FASTA including <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">alternate contigs</a>, corresponding <code>.fai</code> index and <code>.dict</code> dictionary and six BWA-specific index files <code>.alt</code>, <code>.sa</code>, <code>.amb</code>, <code>.bwt</code>, <code>.ann</code> and <code>.pac</code>.</li>
<li>Known variant sites VCFs and corresponding indexes for masking during base quality score recalibration (BQSR). A separate article will document these data resources.</li>
<li>Intervals lists for scattered variant calling with HaplotypeCaller. It is necessary to predefine the calling intervals carefully. In this workflow we use 50 <code>.interval_list</code> files that each contain multiple calling intervals. The calling intervals are an intersection of (i) calling regions of interest and (ii) regions bounded by Ns, otherwise known as gaps in the genome. See the <strong>External Resources</strong> section of <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">Article#7857</a> for an example gap file. Use of these defined intervals has the following benefits. 
<ul>
<li>Avoids calling twice on a locus. This can happen when reads overlapping an interval boundary expand the interval under consideration such that the same variant is called twice. </li>
<li>Makes calls taking into consideration all the available reads that align to the same locus. </li>
</ul></li>
</ul>
<hr />
<p><a name="0"></a></p>
<h1>WORKFLOW definition overview</h1>
<p>Below we see that the workflow name is <strong>PairedEndSingleSampleWorkflow</strong>.</p>
<p>[0.0]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a6/a66257fe80281f12a2e9c6688bd60a.png" />
<p>After the workflow name, the WORKFLOW definition lists the variables that can stand in for files, parameters or even parts of commands within tasks, e.g. the command for BWA alignment (L549). The actual files are given in an accompanying <a href="https://github.com/broadinstitute/wdl/blob/develop/scripts/broad_pipelines/PublicPairedSingleSampleWf_160720.inputs.json"><strong>JSON</strong> file</a>. </p>
<p>[0.1]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/fd/b611570a82f6c18397bc37d0c3471d.png" />
<p>The WORKFLOW definition then outlines the tasks that it will perform. Because tasks may be listed in any order, it is the WORKFLOW definition that defines the order in which steps are run. </p>
<p>Let's break down the workflow into steps and examine their component commands.</p>
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="1"></a></p>
<h2>1. Map with BWA-MEM and merge to create clean BAM</h2>
<h4>This step takes the unaligned BAM, aligns with BWA-MEM, merges information between the unaligned and aligned BAM and fixes tags and sorts the BAM.</h4>
<ul>
<li>1.0: Calls task <a href="#GetBwaVersion">GetBwaVersion</a> to note the version of BWA for use in step [1.3].</li>
<li>1.1:  Defines how to scatter the given unmapped BAM for [1.2], [1.3] and [1.4]. The workflow scatters each <code>unmapped_bam</code> in the list of BAMs given by the variable <code>flowcell_unmapped_bams</code>. The step processes each <code>unmapped_bam</code> in <code>flowcell_unmapped_bams</code> separately in parallel for the three processes. That is, the workflow processes each read group BAM independently for this step.</li>
</ul>
<p>▶︎ Observe the nesting of commands via their relative indentation. Our script writers use these indentations not because they make a difference for Cromwell interpretation but because they allow us human readers to visually comprehend where the scattering applies. In box [1.1] below, we see the scattering defined in L558 applies to processes in boxes [1.2], [1.3] and [1.4] in that the script nests, or indents further in, the commands for these processes within the scattering command. </p>
<ul>
<li>1.2: Calls task <a href="#SamToFastqAndBwaMem">SamToFastqAndBwaMem</a> to map reads to the reference using BWA-MEM. We use <code>bwa_commandline</code> from L549 as the actual command.</li>
<li>1.3: Calls task <a href="#MergeBamAlignment">MergeBamAlignment</a> to merge information between the unaligned and aligned BAM. Both <code>bwa_commandline</code> and <code>bwa_version</code> define elements of the <code>bwamem</code> program group <code>@PG</code> line in the BAM header. The data resulting from this step go on to step [2].</li>
<li>1.4: Calls task <a href="#SortAndFixTags">SortAndFixTags</a> and uses <a href="https://software.broadinstitute.org/wdl/userguide/plumbing#alias">task aliasing</a> to rename it <strong>SortAndFixReadGroupBam</strong>. The consequence of this is that the workflow can then differentiate outputs from those of [2.1]. This task coordinate sorts and indexes the alignments and fixes the <code>NM</code> and <code>UQ</code> tags whose calculations depend on coordinate sort order. This data transformation allows for validation with <a href="https://www.broadinstitute.org/gatk/guide/article?id=7571">ValidateSamFile</a>.</li>
</ul>
<p>[1.0]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/69/e7c6411cec65d5461afa2471cddbff.png" />
<p>[1.1]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/dc/94a73e970d3ffc874eeec647f66ee3.png" />
<p>[1.2]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/83/e617f3821dc523cd55a5caf01b92a8.png" />
<p>[1.3]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/5d/dceade39ab22897b6710b925ad7d10.png" />
<p>[1.4]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/43/c814928512abc21d3cdc057e9bc6ba.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="2"></a></p>
<h2>2. Flag duplicates with MarkDuplicates</h2>
<h4>This step aggregates sample BAMs, flags duplicate sets, fixes tags and coordinate sorts. It starts with the output of [1.3]</h4>
<ul>
<li>2.0: Calls task <a href="#MarkDuplicates">MarkDuplicates</a> to accomplish two aims. First, since MarkDuplicates is given all the files output by the <a href="#MergeBamAlignment">MergeBamAlignment</a> task, and these by design all belong to the same sample, we effectively aggregate the starting BAMs into a single sample-level BAM. Second, MarkDuplicates flags reads it determines are from duplicate inserts with the 0x400 bitwise SAM flag. Because MarkDuplicates sees query-grouped read alignment records from the output of [1.3], it will also mark as duplicate the unmapped mates and supplementary alignments within the duplicate set.  </li>
<li>2.1: Calls task <a href="#SortAndFixTags">SortAndFixTags</a> and renames it as <strong>SortAndFixSampleBam</strong> to differentiate outputs from those of [1.4] that calls the same task. This task coordinate sorts and indexes the alignments and fixes the <code>NM</code> and <code>UQ</code> tags whose calculations depend on coordinate sort order. Resulting data go on to step [3].</li>
</ul>
<p>[2.0]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/7e/5d1fbef5141df46ca39bdd29ebb880.png" />
<p>[2.1]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/72/e482e8034d00ff85c749c1b613e618.png" />
<hr />
<p><a name="3"></a></p>
<h2>3. Base quality score recalibration</h2>
<h4>This step creates intervals for scattering, performs BQSR, merges back the scattered results into a single file and finally compresses the BAM to CRAM format.</h4>
<ul>
<li>3.0: Calls task <a href="#CreateSequenceGroupingTSV">CreateSequenceGroupingTSV</a> to create intervals from the reference genome's <code>.dict</code> dictionary for subsequent use in boxes [3.1] and [3.2]. </li>
<li>3.1: L644 defines scatter intervals as those created in box [3.0] to apply here and for [3.2]. Calls task <a href="#BaseRecalibrator">BaseRecalibrator</a> to act on each interval of the BAM from [2.1] and results in a recalibration table per interval. </li>
<li>3.2: Calls task <a href="#ApplyBQSR">ApplyBQSR</a> to use the <code>GatherBqsrReports.output_bqsr_report</code> from [3.3] and apply the recalibration to the BAM from [2.1] per interval defined by [3.0]. Each resulting recalibrated BAM will contain alignment records from the specified interval including <a href="https://software.broadinstitute.org/gatk/documentation/article?id=6976">unmapped reads from singly mapping pairs</a>. These unmapped records <em>retain</em> SAM alignment information, e.g. mapping contig and coordinate information, but have an asterisk <code>*</code> in the CIGAR field. </li>
<li>3.3: Calls task <a href="#GatherBqsrReports">GatherBqsrReports</a> to consolidate the per interval recalibration tables into a single recalibration table whose sums reflect the consolidated data.</li>
<li>3.4: Calls task <a href="#ApplyBQSR">ApplyBQSR</a> and uses <a href="https://software.broadinstitute.org/wdl/userguide/plumbing#alias">task aliasing</a> to rename it <strong>ApplyBQSRToUnmappedReads</strong>. The consequence of this is that the workflow can then differentiate outputs from those of [3.2]. The step takes as input the BAM from the SortAndFixSampleBam task [2.1] and L697 shows this command runs on unmapped SAM records. These are read pairs the aligner could not map and reads MergeBamAlignment unmapped as contaminants in [1.3] that are at the end of a coordinate-sorted BAM file. The resulting recalibrated BAM contains only such unmapped alignment records. </li>
<li>3.5: Calls task <a href="#GatherBamFiles">GatherBamFiles</a> to concatenate the recalibrated BAMs from [3.2] and [3.4], in order, into a single indexed BAM that retains the header from the first BAM. Resulting data go onto step [4].</li>
<li>3.6: Calls task <a href="#ConvertToCram">ConvertToCram</a> to compress the BAM further in to reference-dependent indexed CRAM format. </li>
</ul>
<p>[3.0]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/b4/5f1dfecbc87ebef57f3406b74ab332.png" />
<p>[3.1]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c4/916e31926d3162d4a945d71b7d04c0.png" />
<p>[3.2]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/26/027f7112aedd910926daa1d1556ae9.png" />
<p>[3.3]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/07/2f96879cde729c2dd6b2f5321a0751.png" />
<p>[3.4]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/a2/5f93440e1c0c49b4afe05ec01e39c5.png" />
<p>[3.5]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/7e/e9aaeaa81ab24dbee1a793a1283d2b.png" />
<p>[3.6]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/39/4a37ba84e9b9f56aeecd8f80b7cd9c.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="4"></a></p>
<h2>4. Call SNP and INDEL variants with HaplotypeCaller</h2>
<h4>This final step uses HaplotypeCaller to call variants over intervals then merges data into a GVCF for the sample, the final output of the workflow.</h4>
<ul>
<li>4.0: Uses scatter intervals defined within the JSON file under <code>scattered_calling_intervals</code> (L728). We use only the <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">primary assembly contigs</a> of GRCh38, grouped into 50 intervals lists, to call variants. Within the GRCh38 intervals lists, the primary assembly's contigs are divided by contiguous regions between regions of Ns. The called task then uses this list of regions to parallelize the task via the <code>-L ${interval_list}</code> option. </li>
</ul>
<p>▶︎ For this pipeline workflow's setup, fifty parallel processes makes sense for a genome of 3 billion basepairs. However, given the same setup, the 50-way split is overkill for a genome of 370 million basepairs as in the case of the <a href="http://www.genomenewsnetwork.org/articles/06_00/puffer_fish.shtml">pufferfish</a>. </p>
<ul>
<li>
<p>4.1: Calls task <a href="#HaplotypeCaller">HaplotypeCaller</a> to call SNP and INDEL variants on the BAM from [3.5] per interval and results in <a href="https://www.broadinstitute.org/gatk/guide/article?id=4017">GVCF format</a> variant calls files per interval. </p>
</li>
<li>
<p>4.2: Calls task <a href="#GatherVCFs">GatherVCFs</a> to merge the per interval GVCFs into one single-sample GVCF. The tool <strong>MergeVcfs</strong> concatenates the GVCFs in the order given by <code>input_vcfs</code> that by this WORKFLOW's design is ordered by contig. </p>
</li>
<li>4.3: Defines files that copy to an output directory if given an OPTIONS JSON file that defines the output directory. If you omit the OPTIONS JSON or omit defining the outputs directory in the OPTIONS JSON, then the workflow skips this step.  </li>
</ul>
<p>[4.0]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/c0/cbee02e274d18689ce0d043678853d.png" />
<p>[4.1]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/8e/ae4111cc2aae71883ff2f960542c00.png" />
<p>[4.2]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ad/8fa051d427ab3dd358380381626567.png" />
<p>[4.3]</p>
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/97/b2d75f340bca3c7ca7a7e9bf9e436f.png" />
<p><a href="#top">back to top</a></p>
<hr />
<p><a name="tasksoverview"></a></p>
<h1>TASK DEFINITIONS</h1>
<h3>GetBwaVersion</h3>
<p>This task obtains the version of BWA  to later notate within the BAM program group (<code>@PG</code>) line.</p>
<p><a name="GetBwaVersion"></a></p>
<pre><code># Get version of BWA
task GetBwaVersion {
  command {
    /usr/gitc/bwa 2&gt;&amp;1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "1 GB"
  }
  output {
    String version = read_string(stdout())
  }
}</code></pre>
<p><a name="SamToFastqAndBwaMem"></a></p>
<h3>SamToFastqAndBwaMem</h3>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ea/4f93dc6b258935c49b1aa0f8a27a01.jpg" height="120"align="left" border="27"/> The input to this task is an unaligned queryname-sorted BAM and the output is an aligned query-grouped BAM. This step pipes three processes: (i) conversion of BAM to FASTQ reads, (ii) [alternate-contig-aware alignment with BWA-MEM and (iii) conversion of SAM to BAM reads. BWA-MEM requires FASTQ reads as input and produces SAM format reads. This task maps the reads using the BWA command defined as a string variable and in this workflow this string is defined in <a href="#0">[0.1]</a>.  </p>
<ul>
<li><em>Dictionary</em> <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7857">Article#7857</a> defines alternate contigs and other reference genome components.</li>
<li><a href="https://software.broadinstitute.org/gatk/guide/article?id=6483#step3D">Step 3D of Tutorial#6483</a> explains the concepts behind piping these processes. </li>
<li>See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=8017">Tutorial#8017</a> for more details on BWA-MEM's alt-aware alignment. </li>
</ul>
<p>The alt-aware alignment depends on use of GRCh38 as the reference, the versions 0.7.13+ of BWA and the presence of BWA's ALT index from <a href="https://github.com/lh3/bwa/tree/master/bwakit">bwa-kit</a>. If the <code>ref_alt</code> ALT index has no content or is not present, then the script exits with an <code>exit 1</code> error. What this means is that this task is only compatible with a reference with ALT contigs and it only runs in an alt-aware manner.</p>
<pre><code># Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File input_bam
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
  File ref_alt

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int disk_size
  Int preemptible_tries

  command &lt;&lt;&lt;
    set -o pipefail
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it,
    if [ -s ${ref_alt} ]; then
      java -Xmx3000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/${bwa_commandline} /dev/stdin -  2&gt; &gt;(tee ${output_bam_basename}.bwa.stderr.log &gt;&amp;2) | \
      samtools view -1 - &gt; ${output_bam_basename}.bam &amp;&amp; \
      grep -m1 "read .* ALT contigs" ${output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      exit 1;
    fi
  &gt;&gt;&gt;
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}</code></pre>
<p><a name="MergeBamAlignment"></a></p>
<h3>MergeBamAlignment</h3>
<p>This step takes an unmapped BAM and the aligned BAM and merges information from each. Reads, sequence and quality information and meta information from the unmapped BAM merge with the alignment information in the aligned BAM. The BWA version the script obtains from task <a href="#GetBwaVersion">GetBwaVersion</a> is used here in the program group (<code>@PG</code>) <code>bwamem</code>. What is imperative for this step, that is implied by the script, is that the sort order of the unmapped and aligned BAMs are identical, i.e. query-group sorted. The BWA-MEM alignment step outputs reads in exactly the same order as they are input and so groups mates, secondary and supplementary alignments together for a given read name. The merging step requires both files maintain this ordering and will produce a final merged BAM in the same query-grouped order given the <code>SORT_ORDER="unsorted"</code> parameter. This has implications for how the <a href="#MarkDuplicates">MarkDuplicates task</a> will flag duplicate sets.</p>
<p>Because the <code>ATTRIBUTES_TO_RETAIN</code> option is set to <code>X0</code>, any aligner-specific tags that are literally <code>X0</code> will carryover to the merged BAM. BWA-MEM does not output such a tag but does output <code>XS</code> and <code>XA</code> tags for suboptimal alignment score and alternative hits, respectively. However, these do not carryover into the merged BAM. Merging retains certain tags from either input BAM (<code>RG</code>, <code>SA</code>, <code>MD</code>, <code>NM</code>, <code>AS</code> and <code>OQ</code> if present), replaces the <code>PG</code> tag as the command below defines and adds new tags (<code>MC</code>, <code>MQ</code> and <code>FT</code>). </p>
<p>▶︎ Note the <code>NM</code> tag values will be incorrect at this point and the <code>UQ</code> tag is absent. Update and addition of these are dependent on coordinate sort order. Specifically, the script uses a separate <a href="#SortAndFixTags">SortAndFixTags</a> task to fix <code>NM</code> tags and add <code>UQ</code> tags. </p>
<p>The <code>UNMAP_CONTAMINANT_READS=true</code> option applies to likely cross-species contamination, e.g. bacterial contamination. MergeBamAlignment identifies reads that are (i) softclipped on both ends and (ii) map with less than 32 basepairs as contaminant. For a similar feature in GATK, see <a href="https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_engine_filters_OverclippedReadFilter.php">OverclippedReadFilter</a>. If MergeBamAlignment determines a read is contaminant, then the mate is also considered contaminant. MergeBamAlignment unmaps the pair of reads by (i) setting the 0x4 flag bit, (ii) replacing column 3's contig name with an asterisk <code>*</code>, (iii) replacing columns 4 and 5 (POS and MAPQ) with zeros, and (iv) adding the <code>FT</code> tag to indicate the reason for unmapping the read, e.g. <code>FT:Z:Cross-species contamination</code>. The records retain their CIGAR strings. Note other processes also use the <code>FT</code> tag, e.g. to indicate reasons for setting the QCFAIL 0x200 flag bit, and will use different tag descriptions. </p>
<pre><code># Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Int disk_size
  Int preemptible_tries

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    java -Xmx3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAP_CONTAMINANT_READS=true
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3500 MB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}</code></pre>
<p><a name="MarkDuplicates"></a></p>
<h3>MarkDuplicates</h3>
<p>This task flags duplicate reads. Because the input is query-group-sorted, MarkDuplicates flags with the 0x400 bitwise SAM flag duplicate primary alignments as well as the duplicate set's secondary and supplementary alignments. Also, for singly mapping mates, duplicate flagging extends to cover unmapped mates. These extensions are features that are only available to query-group-sorted BAMs. </p>
<p>This command uses the <code>ASSUME_SORT_ORDER="queryname"</code> parameter to tell the tool the sort order to expect. Within the context of this workflow, at the point this task is called, we will have avoided any active sorting that would label the BAM header. We know that our original flowcell BAM is queryname-sorted and that BWA-MEM maintains this order to give us query-grouped alignments.</p>
<p>The <code>OPTICAL_DUPLICATE_PIXEL_DISTANCE</code> of 2500 is set for Illumina sequencers that use patterned flowcells to <em>estimate</em> the number of sequencer duplicates. Sequencer duplicates are a subspecies of the duplicates that the tool flags. The Illumina HiSeq X and HiSeq 4000 platforms use patterened flowcells. If <a href="https://software.broadinstitute.org/gatk/guide/article?id=6747#section4">estimating library complexity</a> (see section <em>Duplicate metrics in brief</em>) is important to you, then adjust the <code>OPTICAL_DUPLICATE_PIXEL_DISTANCE</code> appropriately for your sequencer platform. </p>
<p>Finally, in this task and others, we produce an MD5 file with the <code>CREATE_MD5_FILE=true</code> option. This creates a 128-bit hash value using the <a href="https://en.wikipedia.org/wiki/MD5">MD5 algorithm</a> that is to files much like a fingerprint is to an individual. Compare MD5 values to verify data integrity, e.g. after moving or copying large files.</p>
<pre><code># Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Int disk_size

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped, and thus so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"
      CREATE_MD5_FILE=true
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}</code></pre>
<p><a name="SortAndFixTags"></a></p>
<h3>SortAndFixTags</h3>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/ea/4f93dc6b258935c49b1aa0f8a27a01.jpg" height="120"align="left" border="27"/> This task (i) sorts reads by coordinate and then (ii) corrects the NM tag values, adds UQ tags and indexes a BAM. The task pipes the two commands. First, SortSam sorts the records by genomic coordinate using the <code>SORT_ORDER="coordinate"</code> option. Second, SetNmAndUqTags calculates and fixes the UQ and NM tag values in the BAM. Because <code>CREATE_INDEX=true</code>, SetNmAndUqTags creates the <code>.bai</code> index. Again, we create an MD5 file with the <code>CREATE_MD5_FILE=true</code> option.  </p>
<p>As mentioned in the <a href="#MergeBamAlignment">MergeBamAlignment</a> task, tag values dependent on coordinate-sorted records require correction in this separate task given this workflow maintains query-group ordering through the pre-processing steps. </p>
<pre><code># Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
    SortSam \
    INPUT=${input_bam} \
    OUTPUT=/dev/stdout \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false | \
    java -Xmx500m -jar /usr/gitc/picard.jar \
    SetNmAndUqTags \
    INPUT=/dev/stdin \
    OUTPUT=${output_bam_basename}.bam \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true \
    REFERENCE_SEQUENCE=${ref_fasta}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MB"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}</code></pre>
<p><a name="CreateSequenceGroupingTSV"></a></p>
<h3>CreateSequenceGroupingTSV</h3>
<p>This task uses a python script written as a single command using <a href="https://en.wikipedia.org/wiki/Here_document">heredoc syntax</a> to create a list of contig groupings. The workflow uses the intervals to scatter the base quality recalibration step [3] that calls on BaseRecalibrator and ApplyBQSR tasks. </p>
<p>This workflow specifically uses <a href="https://docs.python.org/2/">Python v2.7</a>.  </p>
<p>The input to the task is the reference <code>.dict</code> dictionary that lists contigs. The code takes the information provided by the <code>SN</code> and <code>LN</code> tags of each <code>@SQ</code> line in the dictionary to pair the information in a tuple list. The <code>SN</code> tag names a contig while the <code>LN</code> tag measures the contig length. This list is ordered by descending contig length.  </p>
<p>The contig groupings this command creates is in WDL array format where each line represents a group and each group's members are tab-separated. The command adds contigs to each group from the previously length-sorted list in descending order and caps the sum of member lengths by the first contig's sequence length (the longest contig). This has the effect of somewhat evenly distributing sequence per group. For GRCh38, <code>CreateSequenceGroupingTSV-stdout.log</code> shows 18 such groups. </p>
<p>As the code adds contig names to groups, it adds a <code>:1+</code> to the end of each name. This is to protect the names from downstream tool behavior that removes elements after the last <code>:</code> within a contig name. GRCh38 introduces contig names that include <code>:</code>s and removing the last element make certain contigs indistinguishable from others. With this appendage, we preserve the original contig names through downstream processes. GATK v3.5 and prior versions require this addition. </p>
<pre><code># Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict
  Int preemptible_tries

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.  It outputs to stdout
  # where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command &lt;&lt;&lt;
    python &lt;&lt;CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] &lt;= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    print tsv_string
    CODE
  &gt;&gt;&gt;
  runtime {
    docker: "python:2.7"
    memory: "2 GB"
    preemptible: preemptible_tries
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv(stdout())
  }
}</code></pre>
<p><a name="BaseRecalibrator"></a></p>
<h3>BaseRecalibrator</h3>
<p>The task runs BaseRecalibrator to detect errors made by the sequencer in estimating base quality scores. BaseRecalibrator builds a model of covariation from mismatches in the alignment data while excluding known variant sites and creates a recalibration report for use in the next step. The engine parameter <code>--useOriginalQualities</code> asks BaseRecalibrator to use original sequencer-produced base qualities stored in the <code>OQ</code> tag if present or otherwise use the standard QUAL score. The known sites files should include sites of known common SNPs and INDELs.</p>
<p>This task runs per interval grouping defined by each <code>-L</code> option. The <code>sep</code> in  <code>-L ${sep=" -L " sequence_group_interval}</code> ensures each interval in the _sequence_group<em>interval</em> list is given by the command.</p>
<pre><code># Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx4000m \
      -jar /usr/gitc/GATK4.jar \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${sep=" -knownSites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "6 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}</code></pre>
<p><a name="GatherBqsrReports"></a></p>
<h3>GatherBqsrReports</h3>
<p>This task consolidates the recalibration reports from each sequence group interval into a single report using GatherBqsrReports. </p>
<pre><code># Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx3000m -jar /usr/gitc/GATK4.jar \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}</code></pre>
<p><a name="ApplyBQSR"></a></p>
<h3>ApplyBQSR</h3>
<p>The task uses ApplyBQSR and the recalibration report to correct base quality scores in the BAM. Again, using parallelization, this task applies recalibration for the sequence intervals defined with <code>-L</code>. A resulting recalibrated BAM will contain only reads for the intervals in the applied intervals list.</p>
<pre><code># Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  Int preemptible_tries

  command {
    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
      -jar /usr/gitc/GATK4.jar \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40 \
      --emit_original_quals \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}</code></pre>
<p><a name="GatherBamFiles"></a></p>
<h3>GatherBamFiles</h3>
<p>This task concatenates provided BAMs in order, into a single BAM and retains the header of the first file. For this pipeline, this includes the recalibrated sequence grouped BAMs and the recalibrated unmapped reads BAM. For GRCh38, this makes 19 BAM files that the task concatenates together. The resulting BAM is already in coordinate-sorted order. The task creates a new sequence index and MD5 file for the concatenated BAM. </p>
<pre><code># Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  File input_unmapped_reads_bam
  String output_bam_basename
  Int disk_size
  Int preemptible_tries

  command {
    java -Xmx2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      INPUT=${input_unmapped_reads_bam} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true

    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}</code></pre>
<p><a name="ConvertToCram"></a></p>
<h3>ConvertToCram</h3>
<p>This task compresses a BAM to an even smaller <a href="https://samtools.github.io/hts-specs/CRAMv3.pdf">CRAM format</a> using the <code>-C</code> option of Samtools. The task then indexes the CRAM and renames it from <code>{basename}.cram.crai</code> to <code>{basename}.crai</code>. CRAM is a new format and tools are actively refining features for compatibility. Make sure your tool chain is compatible with CRAM before deleting BAMs. Be aware when using CRAMs that you will have to specify the <em>identical</em> reference genome, not just <em>equivalent</em> reference, with matching MD5 hashes for each contig. These can differ if the capitalization of reference sequences differ.</p>
<pre><code># Convert BAM file to CRAM format
task ConvertToCram {
  File input_bam
  File ref_fasta
  File ref_fasta_index
  String output_basename
  Int disk_size

  # Note that we are not activating pre-emptible instances for this step yet,
  #  but we should if it ends up being fairly quick
  command &lt;&lt;&lt;
      samtools view -C -T ${ref_fasta} ${input_bam} | \
      tee ${output_basename}.cram | \
      md5sum &gt; ${output_basename}.cram.md5 &amp;&amp; \
      samtools index ${output_basename}.cram &amp;&amp; \
      mv ${output_basename}.cram.crai ${output_basename}.crai
  &gt;&gt;&gt;
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "${output_basename}.cram"
    File output_cram_index = "${output_basename}.crai"
    File output_cram_md5 = "››${output_basename}.cram.md5"
  }
}</code></pre>
<p><a name="HaplotypeCaller"></a></p>
<h3>HaplotypeCaller</h3>
<p>This task runs HaplotypeCaller on the recalibrated BAM for given intervals and produces variant calls in <a href="https://www.broadinstitute.org/gatk/guide/article?id=4017">GVCF format</a>. HaplotypeCaller reassembles and realign reads around variants and calls genotypes and genotype likelihoods for single nucleotide polymorphism (SNP) and insertion and deletion (INDELs) variants. Proximal variants are phased. The resulting file is a GZ compressed file, a valid VCF format file with extension <code>.vcf.gz</code>, containing variants for the given interval. </p>
<ul>
<li>The WORKFLOW's <a href="#4">step 4</a> defines any parallelization.</li>
</ul>
<p>The <code>-ERC GVCF</code> or <em>emit reference confidence</em> mode activates two GVCF features. First, for each variant call, we now include a symbolic <code>&lt;NON_REF&gt;</code> <em>non-reference allele</em>. Second, for non-variant regions, we now include <code>&lt;NON_REF&gt;</code> summary blocks as calls. </p>
<ul>
<li>The <code>--max_alternate_alleles</code> is set to three for performance optimization. This does not limit the alleles that are genotyped, only the number of alleles that HaplotypeCaller emits. </li>
<li>Because this WORKFLOW's naming convention does not use the <code>.g.vcf</code> extension, we must specify <code>-variant_index_parameter 128000</code> and <code>-variant_index_type LINEAR</code> to set the correct index strategy for the output GVCF. See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=3893">Article#3893</a> for details.</li>
<li>The command invokes an additional read, the <a href="(https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_engine_filters_OverclippedReadFilter.php)">OverclippedReadFilter</a>, with <code>--read_filter OverclippedRead</code> that removes  reads that are likely from foreign contaminants, e.g. bacterial contamination. The filter define such reads as those that align with less than 30 basepairs and are softclipped on both ends of the read. This option is similar to the <a href="#MergeBamAlignment">MergeBamAlignment task</a>'s <code>UNMAP_CONTAMINANT_READS=true</code> option that unmaps contaminant reads less than 32 basepairs.</li>
</ul>
<pre><code># Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int disk_size
  Int preemptible_tries

  # tried to find lowest memory variable where it would still work, might change once tested on JES
  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${input_bam} \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}</code></pre>
<p><a name="GatherVCFs"></a></p>
<h3>GatherVCFs</h3>
<p>The task uses MergeVcfs to combine multiple VCF files into a single VCF file and index. </p>
<pre><code># Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task GatherVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size
  Int preemptible_tries

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx2g -jar /usr/gitc/picard.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT=${output_vcf_name}
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
}</code></pre>
<p><a href="#top">back to top</a></p>
<hr />