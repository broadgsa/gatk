## Base Quality Score Recalibration (BQSR)

http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr

<p>BQSR stands for Base Quality Score Recalibration. In a nutshell, it is a data pre-processing step that detects systematic errors made by the sequencer when it estimates the quality score of each base call. This document starts with a high-level overview of the purpose of this method; deeper technical are provided further down.</p>
<p>Note that this base recalibration process (BQSR) should not be confused with variant recalibration (VQSR), which is a sophisticated filtering technique applied on the variant callset produced in a later step. <em>The developers who named these methods wish to apologize sincerely to any Spanish-speaking users who might get awfully confused at this point.</em></p>
<hr />
<h3>Wait, what are base quality scores again?</h3>
<p>These scores are per-base estimates of error emitted by the sequencing machines; they express how confident the machine was that it called the correct base each time. For example, let's say the machine reads an A nucleotide, and assigns a quality score of Q20 -- in Phred-scale, that means it's 99% sure it identified the base correctly. This may seem high, but it does mean that we can expect it to be wrong in one case out of 100; so if we have several billion basecalls (we get ~90 billion in a 30x genome), at that rate the machine would make the wrong call in 900 million bases. In practice each basecall gets its own quality score, determined through some dark magic jealously guarded by the manufacturer of the sequencer. </p>
<p>Variant calling algorithms rely heavily on the quality score assigned to the individual base calls in each sequence read. This is because the quality score tells us how much we can trust that particular observation to inform us about the biological truth of the site where that base aligns. If we have a basecall that has a low quality score, that means we're not sure we actually read that A correctly, and it could actually be something else. So we won't trust it as much as other base calls that have higher qualities. In other words we use that score to weigh the evidence that we have for or against a variant allele existing at a particular site. </p>
<h3>Okay, so what is base recalibration?</h3>
<p>Unfortunately the scores produced by the machines are subject to various sources of systematic (non-random) technical error, leading to over- or under-estimated base quality scores in the data. Some of these errors are due to the physics or the chemistry of how the sequencing reaction works, and some are probably due to manufacturing flaws in the equipment.</p>
<p>Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. For example we can identify that, for a given run, whenever we called two A nucleotides in a row, the next base we called had a 1% higher rate of error. So any base call that comes after AA in a read should have its quality score reduced by 1%. We do that over several different covariates (mainly sequence context and position in read, or cycle) in a way that is additive. So the same base may have its quality score increased for one reason and decreased for another.  </p>
<p>This allows us to get more accurate base qualities overall, which in turn improves the accuracy of our variant calls. To be clear, we can't correct the base calls themselves, <em>i.e.</em> we can't determine whether that low-quality A should actually have been a T -- but we can at least tell the variant caller more accurately how far it can trust that A. Note that in some cases we may find that some bases should have a higher quality score, which allows us to rescue observations that otherwise may have been given less consideration than they deserve. Anecdotally my impression is that sequencers are more often over-confident than under-confident, but we do occasionally see runs from sequencers that seemed to suffer from low self-esteem. </p>
<h3>Fantastic! How does it work?</h3>
<p>The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants, then it adjusts the base quality scores in the data based on the model. The known variants are used to mask out bases at sites of real (expected) variation, to avoid counting real variants as errors. Outside of the masked sites, every mismatch is counted as an error. The rest is mostly accounting. </p>
<p>There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes.</p>
<hr />
<h2>More detailed information</h2>
<p>Detailed information about command line options for BaseRecalibrator can be found <a rel="nofollow" class="external text" href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php">here</a>.</p>
<p>The tools in this package recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome.  Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so provides not only more accurate quality scores but also more widely dispersed ones.  The system works on BAM files coming from many sequencing platforms: Illumina, SOLiD, 454, Complete Genomics, Pacific Biosciences, etc. </p>
<p>This process is accomplished by analyzing the covariation among several features of a base. For example: 
</p> 
<ul><li> Reported quality score
</li><li> The position within the read
</li><li> The preceding and current nucleotide (sequencing chemistry effect) observed by the sequencing machine
</li></ul> 
<p>These covariates are then subsequently applied through a piecewise tabular correction to recalibrate the quality scores of all reads in a BAM file. 
</p><p>For example, pre-calibration a file could contain only reported Q25 bases, which seems good.  However, it may be that these bases actually mismatch the reference at a 1 in 100 rate, so are actually Q20.  These higher-than-empirical quality scores provide false confidence in the base calls.  Moreover, as is common with sequencing-by-synthesis machine, base mismatches with the reference occur at the end of the reads more frequently than at the beginning.  Also, mismatches are strongly associated with sequencing context, in that the dinucleotide AC is often much lower quality than TG.  The recalibration tool will not only correct the average Q inaccuracy (shifting from Q25 to Q20) but identify subsets of high-quality bases by separating the low-quality end of read bases AC bases from the high-quality TG bases at the start of the read.  See below for examples of pre and post corrected values.
</p><p>The system was designed for (sophisticated) users to be able to easily add new covariates to the calculations. For users wishing to add their own covariate simply look at QualityScoreCovariate.java for an idea of how to implement the required interface. Each covariate is a Java class which implements the org.broadinstitute.sting.gatk.walkers.recalibration.Covariate interface. Specifically, the class needs to have a getValue method defined which looks at the read and associated sequence context and pulls out the desired information such as machine cycle.
</p> 
<h2><span class="mw-headline" id="Running_the_tools"> Running the tools </span></h2> 
<h3><span class="mw-headline" id="BaseRecalibrator"> BaseRecalibrator </span></h3> 
<p>Detailed information about command line options for BaseRecalibrator can be found <a rel="nofollow" class="external text" href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php">here</a>.
</p><p>This GATK processing step walks over all of the reads in <code>my_reads.bam</code> and tabulates data about the following features of the bases:
</p> 
<ul>
<li>read group the read belongs to</li>
<li>assigned quality score</li>
<li>machine cycle producing this base</li>
<li>current base + previous base (dinucleotide)</li>
</ul> 
<p>For each bin, we count the number of bases within the bin and how often such bases mismatch the reference base, excluding loci known to vary in the population, according to dbSNP.  After running over all reads, BaseRecalibrator produces a file called <code>my_reads.recal_data.grp</code>, which contains the data needed to recalibrate reads.  The format of this GATK report is described below.
</p> 
<h3>Creating a recalibrated BAM</h3> 
<p>To create a recalibrated BAM you can use GATK's PrintReads with the engine on-the-fly recalibration capability. Here is a typical command line to do so:
</p>
<pre> 
java -jar GenomeAnalysisTK.jar \
   -T PrintReads \
   -R reference.fasta \
   -I input.bam \
   -BQSR recalibration_report.grp \
   -o output.bam
</pre> 
<p>After computing covariates in the initial BAM File, we then walk through the BAM file again and rewrite the quality scores (in the QUAL field) using the data in the <code>recalibration_report.grp</code> file, into a new BAM file.  
</p>
<p>This step uses the recalibration table data in recalibration_report.grp produced by BaseRecalibration to recalibrate the quality scores in input.bam, and writing out a new BAM file output.bam with recalibrated QUAL field values.
</p>
<p>Effectively the new quality score is:</p>
<ul>
<li>the sum of the global difference between reported quality scores and the empirical quality</li>
<li>plus the quality bin specific shift</li>
<li>plus the cycle x qual and dinucleotide x qual effect</li>  
</ul>
<p>Following recalibration, the read quality scores are much closer to their empirical scores than before.  This means they can be used in a statistically robust manner for downstream processing, such as SNP calling.  In additional, by accounting for quality changes by cycle and sequence context, we can identify truly high quality bases in the reads, often finding a subset of bases that are Q30 even when no bases were originally labeled as such.
</p> 
<h3>Miscellaneous information</h3> 
<ul><li> The recalibration system is read-group aware.  It separates the covariate data by read group in the recalibration_report.grp file (using @RG tags) and PrintReads will apply this data for each read group in the file.  We routinely process BAM files with multiple read groups.  Please note that the memory requirements scale linearly with the number of read groups in the file, so that files with many read groups could require a significant amount of RAM to store all of the covariate data.
</li>
<li> A critical determinant of the quality of the recalibation is the number of observed bases and mismatches in each bin.  The system will not work well on a small number of aligned reads.  We usually expect well in excess of 100M bases from a next-generation DNA sequencer per read group.  1B bases yields significantly better results.
</li>
<li> Unless your database of variation is so poor and/or variation so common in your organism that most of your mismatches are real snps, you should always perform recalibration on your bam file.  For humans, with dbSNP and now 1000 Genomes available, almost all of the mismatches - even in cancer - will be errors, and an accurate error model (essential for downstream analysis) can be ascertained.
</li>
<li> The recalibrator applies a "yates" correction for low occupancy bins.  Rather than inferring the true Q score from # mismatches / # bases we actually infer it from (# mismatches + 1) / (# bases + 2).  This deals very nicely with overfitting problems, which has only a minor impact on data sets with billions of bases but is critical to avoid overconfidence in rare bins in sparse data.
</li></ul> 
<h2><span class="mw-headline" id="Example_pre_and_post_recalibration_results"> Example pre and post recalibration results </span></h2> 
<ul><li> Recalibration of a lane sequenced at the Broad by an Illumina GA-II in February 2010
</li><li> There is a significant improvement in the accuracy of the base quality scores after applying the GATK recalibration procedure
</li></ul> 
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/d0/d306c3a2d28693598398b8c5443157.png" />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/6f/5309fc58b1e90cfedced982c9cda83.png" />
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/3f/84bd0aa49ae24edf3749dbb6d69cad.png" /> 
<img src="https://us.v-cdn.net/5019796/uploads/FileUpload/16/257337235569ff4ea8f4a05e803b8c.png" />
</p> 
<h2><span class="mw-headline" id="Output"> The output of the BaseRecalibrator </span></h2> 
<ul><li> A Recalibration report containing all the recalibration information for the data
</li></ul> 
<p>Note that the BasRecalibrator no longer produces plots; this is now done by the AnalyzeCovariates tool.</p>
<h3><span class="mw-headline" id="The_Recalibration_Report">The Recalibration Report</span></h3> 
<p>The recalibration report is a [GATKReport](http://gatk.vanillaforums.com/discussion/1244/what-is-a-gatkreport) and not only contains the main result of the analysis, but it is also used as an input to all subsequent analyses on the data. The recalibration report contains the following 5 tables:
</p> 
<ul><li> Arguments Table -- a table with all the arguments and its values 
</li><li> Quantization Table
</li><li> ReadGroup Table
</li><li> Quality Score Table
</li><li> Covariates Table
</li></ul> 
<h4><span class="mw-headline" id="Arguments_Table">Arguments Table</span></h4> 
<p>This is the table that contains all the arguments used to run BQSRv2 for this dataset. This is important for the on-the-fly recalibration step to use the same parameters used in the recalibration step (context sizes, covariates, ...).
</p>
<p>Example Arguments table:</p> 
<pre> 
#:GATKTable:true:1:17::;
#:GATKTable:Arguments:Recalibration argument collection values used in this run
Argument                    Value
covariate                   null
default_platform            null
deletions_context_size      6
force_platform              null
insertions_context_size     6
...
</pre> 
<h4><span class="mw-headline" id="Quantization_Table">Quantization Table</span></h4> 
<p>The GATK offers native support to quantize base qualities. The GATK quantization procedure uses a statistical approach to determine the best binning system that minimizes the error introduced by amalgamating the different qualities present in the specific dataset. When running BQSRv2, a table with the base counts for each base quality is generated and a 'default' quantization table is generated. This table is a required parameter for any other tool in the GATK if you want to quantize your quality scores.
</p>
<p>The default behavior (currently) is to use no quantization when performing on-the-fly recalibration. You can override this by using the engine argument -qq. With -qq 0 you don't quantize qualities, or -qq N you recalculate the quantization bins using N bins on the fly.  Note that quantization is completely experimental now and we do not recommend using it unless you are a super advanced user.
</p>
<p>Example Arguments table:
</p> 
<pre> 
#:GATKTable:true:2:94:::;
#:GATKTable:Quantized:Quality quantization map
QualityScore  Count        QuantizedScore
0                     252               0
1                   15972               1
2                  553525               2
3                 2190142               9
4                 5369681               9
9                83645762               9
...
</pre> 
<h4><span class="mw-headline" id="ReadGroup_Table">ReadGroup Table</span></h4> 
<p>This table contains the empirical quality scores for each read group, for mismatches insertions and deletions. This is not different from the table used in the old table recalibration walker.
</p> 
<pre> 
#:GATKTable:false:6:18:%s:%s:%.4f:%.4f:%d:%d:;
#:GATKTable:RecalTable0:
ReadGroup  EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors
SRR032768  D                   40.7476             45.0000    2642683174    222475
SRR032766  D                   40.9072             45.0000    2630282426    213441
SRR032764  D                   40.5931             45.0000    2919572148    254687
SRR032769  D                   40.7448             45.0000    2850110574    240094
SRR032767  D                   40.6820             45.0000    2820040026    241020
SRR032765  D                   40.9034             45.0000    2441035052    198258
SRR032766  M                   23.2573             23.7733    2630282426  12424434
SRR032768  M                   23.0281             23.5366    2642683174  13159514
SRR032769  M                   23.2608             23.6920    2850110574  13451898
SRR032764  M                   23.2302             23.6039    2919572148  13877177
SRR032765  M                   23.0271             23.5527    2441035052  12158144
SRR032767  M                   23.1195             23.5852    2820040026  13750197
SRR032766  I                   41.7198             45.0000    2630282426    177017
SRR032768  I                   41.5682             45.0000    2642683174    184172
SRR032769  I                   41.5828             45.0000    2850110574    197959
SRR032764  I                   41.2958             45.0000    2919572148    216637
SRR032765  I                   41.5546             45.0000    2441035052    170651
SRR032767  I                   41.5192             45.0000    2820040026    198762
</pre> 
<h4><span class="mw-headline" id="Quality_Score_Table">Quality Score Table</span></h4> 
<p>This table contains the empirical quality scores for each read group and original quality score, for mismatches insertions and deletions. This is not different from the table used in the old table recalibration walker.
</p> 
<pre> 
#:GATKTable:false:6:274:%s:%s:%s:%.4f:%d:%d:;
#:GATKTable:RecalTable1:
ReadGroup  QualityScore  EventType  EmpiricalQuality  Observations  Errors
SRR032767            49  M                   33.7794          9549        3
SRR032769            49  M                   36.9975          5008        0
SRR032764            49  M                   39.2490          8411        0
SRR032766            18  M                   17.7397      16330200   274803
SRR032768            18  M                   17.7922      17707920   294405
SRR032764            45  I                   41.2958    2919572148   216637
SRR032765             6  M                    6.0600       3401801   842765
SRR032769            45  I                   41.5828    2850110574   197959
SRR032764             6  M                    6.0751       4220451  1041946
SRR032767            45  I                   41.5192    2820040026   198762
SRR032769             6  M                    6.3481       5045533  1169748
SRR032768            16  M                   15.7681      12427549   329283
SRR032766            16  M                   15.8173      11799056   309110
SRR032764            16  M                   15.9033      13017244   334343
SRR032769            16  M                   15.8042      13817386   363078
...
</pre> 
<h4><span class="mw-headline" id="Covariates_Table">Covariates Table</span></h4> 
<p>This table has the empirical qualities for each covariate used in the dataset. The default covariates are cycle and context. In the current implementation, context is of a fixed size (default 6). Each context and each cycle will have an entry on this table stratified by read group and original quality score.
</p> 
<pre> 
#:GATKTable:false:8:1003738:%s:%s:%s:%s:%s:%.4f:%d:%d:;
#:GATKTable:RecalTable2:
ReadGroup  QualityScore  CovariateValue  CovariateName  EventType  EmpiricalQuality  Observations  Errors
SRR032767            16  TACGGA          Context        M                   14.2139           817      30
SRR032766            16  AACGGA          Context        M                   14.9938          1420      44
SRR032765            16  TACGGA          Context        M                   15.5145           711      19
SRR032768            16  AACGGA          Context        M                   15.0133          1585      49
SRR032764            16  TACGGA          Context        M                   14.5393           710      24
SRR032766            16  GACGGA          Context        M                   17.9746          1379      21
SRR032768            45  CACCTC          Context        I                   40.7907        575849      47
SRR032764            45  TACCTC          Context        I                   43.8286        507088      20
SRR032769            45  TACGGC          Context        D                   38.7536         37525       4
SRR032768            45  GACCTC          Context        I                   46.0724        445275      10
SRR032766            45  CACCTC          Context        I                   41.0696        575664      44
SRR032769            45  TACCTC          Context        I                   43.4821        490491      21
SRR032766            45  CACGGC          Context        D                   45.1471         65424       1
SRR032768            45  GACGGC          Context        D                   45.3980         34657       0
SRR032767            45  TACGGC          Context        D                   42.7663         37814       1
SRR032767            16  AACGGA          Context        M                   15.9371          1647      41
SRR032764            16  GACGGA          Context        M                   18.2642          1273      18
SRR032769            16  CACGGA          Context        M                   13.0801          1442      70
SRR032765            16  GACGGA          Context        M                   15.9934          1271      31
...
</pre> 
<h2><span class="mw-headline" id="Troubleshooting"> Troubleshooting </span></h2> 
<p><strong>The memory requirements of the recalibrator will vary based on the type of JVM running the application and the number of read groups in the input bam file.</strong></p>
<p>If the application reports 'java.lang.OutOfMemoryError: Java heap space', increase the max heap size provided to the JVM by adding ' -Xmx????m' to the jvm_args variable in RecalQual.py, where '????' is the maximum available memory on the processing computer.</p>
<p><strong>I've tried recalibrating my data using a downloaded file, such as NA12878 on 454, and apply the table to any of the chromosome BAM files always fails due to hitting my memory limit. I've tried giving it as much as 15GB but that still isn't enough.</strong></p>
<p>All of our big merged files for 454 are running with -Xmx16000m arguments to the JVM -- it's enough to process all of the files.  32GB might make the 454 runs a lot faster though.</p>
<p><strong>I have a recalibration file calculated over the entire genome (such as for the 1000 genomes trio) but I split my file into pieces (such as by chromosome).  Can the recalibration tables safely be applied to the per chromosome BAM files?</strong></p>
<p>Yes they can.  The original tables needed to be calculated over the whole genome but they can be applied to each piece of the data set independently.</p>
<p><strong>I'm working on a genome that doesn't really have a good SNP database yet. I'm wondering if it still makes sense to run base quality score recalibration without known SNPs.</strong></p>
<p>The base quality score recalibrator treats every reference mismatch as indicative of machine error. True polymorphisms are legitimate mismatches to the reference and shouldn't be counted against the quality of a base. We use a database of known polymorphisms to skip over most polymorphic sites. Unfortunately without this information the data becomes almost completely unusable since the quality of the bases will be inferred to be much much lower than it actually is as a result of the reference-mismatching SNP sites.</p>
<p>However, all is not lost if you are willing to experiment a bit. You can bootstrap a database of known SNPs. Here's how it works: </p>
<ul>
<li>First do an initial round of SNP calling on your original, unrecalibrated data. </li>
<li>Then take the SNPs that you have the highest confidence in and use that set as the database of known SNPs by feeding it as a VCF file to the base quality score recalibrator.</li>
<li>Finally, do a real round of SNP calling with the recalibrated data. These steps could be repeated several times until convergence.</li>
</ul>
<h3><span class="mw-headline" id="Downsampling_to_reduce_run_time"> Downsampling to reduce run time </span></h3> 
<p>For users concerned about run time please note this small analysis below showing the approximate number of reads per read group that are required to achieve a given level of recalibration performance. The analysis was performed with 51 base pair Illumina reads on pilot data from the 1000 Genomes Project. Downsampling can be achieved by specifying a genome interval using the -L option. For users concerned only with recalibration accuracy please disregard this plot and continue to use all available data when generating the recalibration table.
</p>
<p><img src="https://us.v-cdn.net/5019796/uploads/FileUpload/bf/48a038896124ddc734d09dabeb7cd4.png" />
</p> 