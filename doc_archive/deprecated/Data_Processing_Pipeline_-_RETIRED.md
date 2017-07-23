## Data Processing Pipeline - RETIRED

http://gatkforums.broadinstitute.org/gatk/discussion/41/data-processing-pipeline-retired

<h3>Please note that the DataProcessingPipeline qscript is no longer available. We are looking into the possibility of producing some new Qscripts that will be more appropriate for sharing with the public.</h3>
<p><em>The DPP script was only provided has an example, but many people were using it &quot;out of the box&quot; without properly understanding how it works. In order to protect users from mishandling this tool, and to decrease our support burden, we have taken the difficult decision of removing the script from our public repository. If you would like to put together your own version of the DPP, please have a look at our other example scripts to understand how Qscripts work, and read the Best Practices documentation to understand what are the processing steps and what parameters you need to set/adjust.</em></p>
<h2>Data Processing Pipeline</h2>
<p>The Data Processing Pipeline is a Queue script designed to take BAM files from the NGS machines to <em>analysis ready</em> BAMs for the GATK. </p>
<h3>Introduction</h3>
<p>Reads come off the sequencers in a raw state that is not suitable for analysis using the GATK. In order to prepare the dataset, one must perform the steps described <a href="http://www.broadinstitute.org/gatk/guide/topic?name=best-practices">here</a>. This pipeline performs the following steps: indel cleaning, duplicate marking and base score recalibration, following the GSA's latest definition of best practices. The product of this pipeline is a set of <em>analysis ready</em> BAM files (one per sample sequenced).</p>
<h3>Requirements</h3>
<p>This pipeline is a <a href="http://www.broadinstitute.org/gatk/guide/article?id=1306">Queue</a> script that uses tools from the GATK, <a href="http://picard.sourceforge.net/">Picard</a> and <a href="http://bio-bwa.sourceforge.net/">BWA</a> (optional) software suites which are all freely available through their respective websites. Queue is a GATK companion that is included in the GATK package.</p>
<p><strong>Warning:</strong> This pipeline was designed specifically to handle the Broad Institute's main sequencing pipeline with Illumina BAM files and BWA alignment. The GSA cannot support its use for other types of datasets. It is possible however, with some effort, to modify it for your needs.</p>
<h3>Command-line arguments</h3>
<h4>Required Parameters</h4>
<table border="1" cellpadding="2" width="100%">
<tr>
<th scope="col" width="15%"> Argument (short-name)
</th>
<th scope="col" width="25%"> Argument (long-name)
</th>
<th scope="col"> Description
</th></tr>
<tr>
<td> -i &lt;BAM file / BAM list&gt; </td>
<td> --input &lt;BAM file / BAM list&gt; </td>
<td> input BAM file - or list of BAM files.
</td></tr>
<tr>
<td> -R &lt;fasta&gt; </td>
<td> --reference &lt;fasta&gt; </td>
<td> Reference fasta file.
</td></tr>
<tr>
<td> -D &lt;vcf&gt; </td>
<td> --dbsnp &lt;dbsnp vcf&gt; </td>
<td> dbsnp ROD to use (must be in VCF format).
</td></tr></table>
<h4>Optional Parameters</h4>
<table border="1" cellpadding="2" width="100%">
<tr>
<th scope="col" width="15%"> Argument (short-name)
</th>
<th scope="col" width="25%"> Argument (long-name)
</th>
<th scope="col"> Description
</th></tr>
<tr>
<td> -indels &lt;vcf&gt; </td>
<td> --extra_indels &lt;vcf&gt; </td>
<td> VCF files to use as reference indels for Indel Realignment.
</td></tr>
<tr>
<td> -bwa &lt;path&gt; </td>
<td> --path_to_bwa &lt;path&gt; </td>
<td> The path to the binary of bwa (usually BAM files have already been mapped - but if you want to remap this is the option)
</td></tr>
<tr>
<td> -outputDir &lt;path&gt; </td>
<td> --output_directory &lt;path&gt; </td>
<td> Output path for the processed BAM files.
</td></tr>
<tr>
<td> -L &lt;GATK interval string&gt; </td>
<td> --gatk_interval_string &lt;GATK interval string&gt; </td>
<td> the -L interval string to be used by GATK - output bams at interval only
</td></tr>
<tr>
<td> -intervals &lt;GATK interval file&gt; </td>
<td> --gatk_interval_file &lt;GATK interval file&gt; </td>
<td> an <i>intervals</i> file to be used by GATK - output bams at intervals
</td></tr></table>
<h4>Modes of Operation (also optional parameters)</h4>
<table border="1" cellpadding="2" width="100%">
<tr>
<th scope="col" width="15%"> Argument (short-name)
</th>
<th scope="col" width="25%"> Argument (long-name)
</th>
<th scope="col"> Description
</th></tr>
<tr>
<td> -p &lt;name&gt; </td>
<td> --project &lt;name&gt; </td>
<td> the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam
</td></tr>
<tr>
<td> -knowns </td>
<td> --knowns_only </td>
<td> Perform cleaning on knowns only.
</td></tr>
<tr>
<td> -sw </td>
<td> --use_smith_waterman </td>
<td> Perform cleaning using Smith Waterman
</td></tr>
<tr>
<td> -bwase </td>
<td> --use_bwa_single_ended </td>
<td> Decompose input BAM file and fully realign it using BWA and assume Single Ended reads
</td></tr>
<tr>
<td> -bwape </td>
<td> --use_bwa_pair_ended </td>
<td> Decompose input BAM file and fully realign it using BWA and assume Pair Ended reads
</td></tr></table>
<h2>The Pipeline</h2>
<p>Data processing pipeline of the best practices for raw data processing, from sequencer data (fastq files) to analysis read reads (bam file):</p>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/55/0a67f9e1b7962a14c422e993f34643.jpeg" alt="the data processing pipeline" /></p>
<p>Following the group's Best Practices definition, the data processing pipeline does all the processing at the sample level. There are two high-level parts of the pipeline:</p>
<h3>BWA alignment</h3>
<p>This option is for datasets that have already been processed using a different pipeline or different criteria, and you want to reprocess it using this pipeline. One example is a BAM file that has been processed at the lane level, or did not perform some of the best practices steps of the current pipeline. By using the optional BWA stage of the processing pipeline, your BAM file will be realigned from scratch before creating sample level bams and entering the pipeline.</p>
<h3>Sample Level Processing</h3>
<p>This is the where the pipeline applies its main procedures: Indel Realignment and Base Quality Score Recalibration. </p>
<h4>Indel Realignment</h4>
<p>This is a two step process. First we create targets using the Realigner Target Creator (either for knowns only, or including data indels), then we realign the targets using the Indel Realigner (see [Local realignment around indels]) with an optional smith waterman realignment. The Indel Realigner also fixes mate pair information for reads that get realigned.</p>
<h4>Base Quality Score Recalibration</h4>
<p>This is a crucial step that re-adjusts the quality score using statistics based on several different covariates. In this pipeline we utilize four: Read Group Covariate, Quality Score Covariate, Cycle Covariate, Dinucleotide Covariate</p>
<h3>The Outputs</h3>
<p>The Data Processing Pipeline produces 3 types of output for each file: a fully processed bam file, a validation report on the input bam and output bam files, a analysis before and after base quality score recalibration. If you look at the pipeline flowchart, the grey boxes indicate processes that generate an output. </p>
<h4>Processed Bam File</h4>
<p>The final product of the pipeline is one BAM file per sample in the dataset. It also provides one BAM list with all the bams in the dataset. This file is named &lt;project name&gt;.cohort.list, and each sample bam file has the name &lt;project name&gt;.&lt;sample name&gt;.bam. The sample names are extracted from the input BAM headers, and the project name is provided as a parameter to the pipeline.</p>
<h4>Validation Files</h4>
<p>We validate each unprocessed sample level BAM file and each final processed sample level BAM file. The validation is performed using <a href="http://picard.sourceforge.net/">Picard</a>'s ValidateSamFile. Because the parameters of this validation are very strict, we don't enforce that the input BAM has to pass all validation, but we provide the log of the validation as an informative companion to your input. The validation file is named&#160;: &lt;project name&gt;.&lt;sample name&gt;.pre.validation and &lt;project name&gt;.&lt;sample name&gt;.post.validation.</p>
<p>Notice that even if your BAM file fails validation, the pipeline can still go through successfully. The validation is a strict report on how your BAM file is looking. Some errors are not critical, but the output files (both pre.validation and post.validation) should give you some input on how to make your dataset better organized in the BAM format.</p>
<h4>Base Quality Score Recalibration Analysis</h4>
<p>PDF plots of the base qualities are generated before and after recalibration for further analysis on the impact of recalibrating the base quality scores in each sample file. These graphs are explained in detail <a href="http://www.broadinstitute.org/gatk/guide/article?id=44">here</a>. The plots are created in directories named&#160;: &lt;project name&gt;.&lt;sample name&gt;.pre and &lt;project name&gt;.&lt;sample name&gt;.post.</p>
<h3>Examples</h3>
<ol>
<li>
<p>Example script that runs the data processing pipeline with its standard parameters and uses LSF for scatter/gathering (without bwa)</p>
<p>java \
-Xmx4g \
-Djava.io.tmpdir=/path/to/tmpdir \
-jar path/to/GATK/Queue.jar \
-S path/to/DataProcessingPipeline.scala \
-p myFancyProjectName \
-i myDataSet.list \
-R reference.fasta \
-D dbSNP.vcf \
-run</p>
</li>
<li>
<p>Performing realignment and the full data processing pipeline in one pair-ended bam file</p>
<p>java \
-Xmx4g \
-Djava.io.tmpdir=/path/to/tmpdir \
-jar path/to/Queue.jar \
-S path/to/DataProcessingPipeline.scala \
-bwa path/to/bwa \
-i test.bam \
-R reference.fasta \
-D dbSNP.vcf \
-p myProjectWithRealignment \
-bwape \
-run</p>
</li>
</ol>