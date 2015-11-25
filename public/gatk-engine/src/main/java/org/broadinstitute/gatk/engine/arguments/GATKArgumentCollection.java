/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.arguments;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport;
import org.broadinstitute.gatk.engine.samples.PedigreeValidationType;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.baq.BAQ;
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType;
import org.broadinstitute.gatk.engine.GATKVCFUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * @author aaron
 * @version 1.0
 */
public class GATKArgumentCollection {

    // the default value of the stop of the expanded window
    public static final int DEFAULT_REFERENCE_WINDOW_STOP = 0;

    /** the constructor */
    public GATKArgumentCollection() {
    }

    // parameters and their defaults
    /**
     * An input file containing sequence data mapped to a reference, in BAM or CRAM format, or a text file containing a
     * list of input files (with extension .list). Note that the GATK requires an accompanying .bai index for each BAM
     * or CRAM file. Please see our online documentation for more details on input formatting requirements.
     */
    @Input(fullName = "input_file", shortName = "I", doc = "Input file containing sequence data (BAM or CRAM)", required = false)
    public List<String> samFiles = new ArrayList<>();

    @Advanced
    @Argument(fullName = "showFullBamList",doc="Emit a log entry (level INFO) containing the full list of sequence data files to be included in the analysis (including files inside .bam.list or .cram.list files).")
    public Boolean showFullBamList = false;

    @Advanced
    @Argument(fullName = "read_buffer_size", shortName = "rbs", doc="Number of reads per SAM file to buffer in memory", required = false, minValue = 0)
    public Integer readBufferSize = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // GATKRunReport options
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * By default, GATK generates a run report that is uploaded to a cloud-based service. This report contains basic
     * statistics about the run (which tool was used, whether the run was successful etc.) that help us for debugging
     * and development. Up to version 3.3-0 the run report contains a record of the username and hostname associated
     * with the run, but it does **NOT** contain any information that could be used to identify patient data.
     * Nevertheless, if your data is subject to stringent confidentiality clauses (no outside communication) or if your
     * run environment is not connected to the internet, you can disable the reporting system by seeting this option to
     * "NO_ET". You will also need to request a key using the online request form on our website (see FAQs).
     */
    @Argument(fullName = "phone_home", shortName = "et", doc="Run reporting mode", required = false)
    public GATKRunReport.PhoneHomeOption phoneHomeType = GATKRunReport.PhoneHomeOption.AWS;
    /**
     * Please see the "phone_home" argument below and the online documentation FAQs for more details on the key system
     * and how to request a key.
     */
    @Argument(fullName = "gatk_key", shortName = "K", doc="GATK key file required to run with -et NO_ET", required = false)
    public File gatkKeyFile = null;

    /**
     * The GATKRunReport supports tagging GATK runs with an arbitrary tag that can be
     * used to group together runs during later analysis (as of GATK 2.2) .  One use of this capability is to tag
     * runs as GATK performance tests, so that the performance of the GATK over time can be assessed from the logs
     * directly.
     *
     * Note that the tags do not conform to any ontology, so you are free to use any tags that you might find
     * meaningful.
     */
    @Argument(fullName = "tag", shortName = "tag", doc="Tag to identify this GATK run as part of a group of runs", required = false)
    public String tag = "NA";

    // --------------------------------------------------------------------------------------------------------------
    //
    // General features
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Reads that fail the specified filters will not be used in the analysis. Multiple filters can be specified separately,
     * e.g. you can do -rf MalformedRead -rf BadCigar and so on. Available read filters are listed in the online tool
     * documentation. Note that the read name format is e.g. MalformedReadFilter, but at the command line the filter
     * name should be given without the Filter suffix; e.g. -rf MalformedRead (NOT -rf MalformedReadFilter, which is not
     * recognized by the program). Note also that some read filters are applied by default for some analysis tools; this
     * is specified in each tool's documentation. The default filters can only be disabled if they are DisableableReadFilters.
     */
    @Argument(fullName = "read_filter", shortName = "rf", doc = "Filters to apply to reads before analysis", required = false)
    public final List<String> readFilters = new ArrayList<>();

    @Argument(fullName = "disable_read_filter", shortName = "drf", doc = "Read filters to disable", required = false)
    public final List<String> disabledReadFilters = new ArrayList<>();

    @ArgumentCollection
    public IntervalArgumentCollection intervalArguments = new IntervalArgumentCollection();
    /**
     * The reference genome against which the sequence data was mapped. The GATK requires an index file and a dictionary
     * file accompanying the reference (please see the online documentation FAQs for more details on these files). Although
     * this argument is indicated as being optional, almost all GATK tools require a reference in order to run.
     * Note also that while GATK can in theory process genomes from any organism with any number of chromosomes or contigs,
     * it is not designed to process draft genome assemblies and performance will decrease as the number of contigs in
     * the reference increases. We strongly discourage the use of unfinished genome assemblies containing more than a few
     * hundred contigs. Contig numbers in the thousands will most probably cause memory-related crashes.
     */
    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;
    /**
     * If this flag is enabled, the random numbers generated will be different in every run, causing GATK to behave non-deterministically.
     */
    @Argument(fullName = "nonDeterministicRandomSeed", shortName = "ndrs", doc = "Use a non-deterministic random seed", required = false)
    public boolean nonDeterministicRandomSeed = false;
    /**
     * To be used in the testing framework where dynamic parallelism can result in differing numbers of calls to the random generator.
     */
    @Hidden
    @Argument(fullName = "disableDithering",doc="Completely eliminates randomized dithering from rank sum tests.")
    public boolean disableDithering = false;
    /**
     * This will truncate the run but without exiting with a failure. By default the value is interpreted in minutes, but this can be changed with the maxRuntimeUnits argument.
     */
    @Argument(fullName = "maxRuntime", shortName = "maxRuntime", doc="Stop execution cleanly as soon as maxRuntime has been reached", required = false)
    public long maxRuntime = GenomeAnalysisEngine.NO_RUNTIME_LIMIT;

    @Argument(fullName = "maxRuntimeUnits", shortName = "maxRuntimeUnits", doc="Unit of time used by maxRuntime", required = false)
    public TimeUnit maxRuntimeUnits = TimeUnit.MINUTES;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Downsampling Arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * There are several ways to downsample reads, i.e. to remove reads from the pile of reads that will be used for analysis.
     * See the documentation of the individual downsampling options for details on how they work. Note that many GATK tools
     * specify a default downsampling type and target, but this behavior can be overridden from the command line using the
     * downsampling arguments.
     */
    @Argument(fullName = "downsampling_type", shortName="dt", doc="Type of read downsampling to employ at a given locus", required = false)
    public DownsampleType downsamplingType = null;
    /**
     * Reads will be downsampled so the specified fraction remains; e.g. if you specify -dfrac 0.25, three-quarters of
     * the reads will be removed, and the remaining one quarter will be used in the analysis. This method of downsampling
     * is truly unbiased and random. It is typically used to simulate the effect of generating different amounts of
     * sequence data for a given sample. For example, you can use this in a pilot experiment to evaluate how much target
     * coverage you need to aim for in order to obtain enough coverage in all loci of interest.
     */
    @Argument(fullName = "downsample_to_fraction", shortName = "dfrac", doc = "Fraction of reads to downsample to", required = false, minValue = 0.0, maxValue = 1.0)
    public Double downsampleFraction = null;

    /**
     * The principle of this downsampling type is to downsample reads to a given capping threshold coverage. Its purpose is to
     * get rid of excessive coverage, because above a certain depth, having additional data is not informative and imposes
     * unreasonable computational costs. The downsampling process takes two different forms depending on the type of
     * analysis it is used with.
     *
     * For locus-based traversals (LocusWalkers like UnifiedGenotyper and ActiveRegionWalkers like HaplotypeCaller),
     * downsample_to_coverage controls the maximum depth of coverage at each locus. For read-based traversals
     * (ReadWalkers like BaseRecalibrator), it controls the maximum number of reads sharing the same alignment start
     * position. For ReadWalkers you will typically need to use much lower dcov values than you would with LocusWalkers
     * to see an effect. Note that this downsampling option does not produce an unbiased random sampling from all available
     * reads at each locus: instead, the primary goal of the to-coverage downsampler is to maintain an even representation
     * of reads from all alignment start positions when removing excess coverage. For a truly unbiased random sampling of
     * reads, use -dfrac instead. Also note that the coverage target is an approximate goal that is not guaranteed to be
     * met exactly: the downsampling algorithm will under some circumstances retain slightly more or less coverage than
     * requested.
     */
    @Argument(fullName = "downsample_to_coverage", shortName = "dcov",
              doc = "Target coverage threshold for downsampling to coverage",
              required = false, minValue = 0)
    public Integer downsampleCoverage = null;

    /**
     * Gets the downsampling method explicitly specified by the user. If the user didn't specify
     * a default downsampling mechanism, return the default.
     * @return The explicitly specified downsampling mechanism, or the default if none exists.
     */
    public DownsamplingMethod getDownsamplingMethod() {
        if ( downsamplingType == null && downsampleFraction == null && downsampleCoverage == null )
            return null;

        return new DownsamplingMethod(downsamplingType, downsampleCoverage, downsampleFraction);
    }

    /**
     * Set the downsampling method stored in the argument collection so that it is read back out when interrogating the command line arguments.
     * @param method The downsampling mechanism.
     */
    public void setDownsamplingMethod(DownsamplingMethod method) {
        if (method == null)
            throw new IllegalArgumentException("method is null");

        downsamplingType = method.type;
        downsampleCoverage = method.toCoverage;
        downsampleFraction = method.toFraction;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // BAQ arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    @Advanced
    @Argument(fullName = "baq", shortName="baq", doc="Type of BAQ calculation to apply in the engine", required = false)
    public BAQ.CalculationMode BAQMode = BAQ.CalculationMode.OFF;
    /**
     *  Phred-scaled gap open penalty for BAQ calculation. Although the default value is 40, a value of 30 may be better for whole genome call sets.
     */
    @Advanced
    @Argument(fullName = "baqGapOpenPenalty", shortName="baqGOP", doc="BAQ gap open penalty", required = false, minValue = 0)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    // --------------------------------------------------------------------------------------------------------------
    //
    // refactor NDN cigar string arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Some RNAseq aligners that use a known transcriptome resource (such as TopHat2) produce NDN elements in read CIGARS
     * when a small exon is entirely deleted during transcription, which ends up looking like [exon1]NDN[exon3]. These
     * rarely happen, but when they do they cause GATK to fail with an error. Setting this flag tells the GATK to
     * reduce "NDN" to a simpler CIGAR representation with one N element (with total length of the three refactored
     * elements). From the point of view of variant calling, there is no meaningful difference between the two
     * representations.
     */
    @Argument(fullName = "refactor_NDN_cigar_string", shortName = "fixNDN", doc = "Reduce NDN elements in CIGAR string", required = false)
    public boolean REFACTOR_NDN_CIGAR_READS = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // quality encoding checking arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * By default the GATK assumes that base quality scores start at Q0 == ASCII 33 according to the SAM specification.
     * However, encoding in some datasets (especially older Illumina ones) starts at Q64. This argument will fix the
     * encodings on the fly (as the data is read in) by subtracting 31 from every quality score. Note that this argument should
     * NEVER be used by default; you should only use it when you have confirmed that the quality scores in your data are
     * not in the correct encoding.
     */
    @Argument(fullName = "fix_misencoded_quality_scores", shortName="fixMisencodedQuals", doc="Fix mis-encoded base quality scores", required = false)
    public boolean FIX_MISENCODED_QUALS = false;
    /**
     * This flag tells GATK to ignore warnings when encountering base qualities that are too high and that seemingly
     * indicate a problem with the base quality encoding of the BAM or CRAM file. You should only use this if you really
     * know what you are doing; otherwise you could seriously mess up your data and ruin your analysis.
     */
    @Argument(fullName = "allow_potentially_misencoded_quality_scores", shortName="allowPotentiallyMisencodedQuals", doc="Ignore warnings about base quality score encoding", required = false)
    public boolean ALLOW_POTENTIALLY_MISENCODED_QUALS = false;
    /**
     * This flag tells GATK to use the original base qualities (that were in the data before BQSR/recalibration) which
     * are stored in the OQ tag, if they are present, rather than use the post-recalibration quality scores. If no OQ
     * tag is present for a read, the standard qual score will be used.
     */
    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "Use the base quality scores from the OQ tag", required=false)
    public Boolean useOriginalBaseQualities = false;
    /**
     * If reads are missing some or all base quality scores, this value will be used for all base quality scores.
     * By default this is set to -1 to disable default base quality assignment.
     */
    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "Assign a default base quality", required=false, minValue = 0, maxValue = Byte.MAX_VALUE)
    public byte defaultBaseQualities = -1;

    // --------------------------------------------------------------------------------------------------------------
    //
    // performance log arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * The file name for the GATK performance log output, or null if you don't want to generate the
     * detailed performance logging table.  This table is suitable for importing into R or any
     * other analysis software that can read tsv files.
     */
    @Argument(fullName = "performanceLog", shortName="PF", doc="Write GATK runtime performance log to this file", required = false)
    public File performanceLog = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // BQSR arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Enables on-the-fly recalibrate of base qualities, intended primarily for use with BaseRecalibrator and PrintReads
     * (see Best Practices workflow documentation). The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input BAM(s)
     * or CRAM(s).
     */
    @Input(fullName="BQSR", shortName="BQSR", required=false, doc="Input covariates table file for on-the-fly base quality score recalibration")
    public File BQSR_RECAL_FILE = null;

    /**
     * Turns on the base quantization module. It requires a recalibration report (-BQSR).
     *
     * A value of 0 here means "do not quantize".
     * Any value greater than zero will be used to recalculate the quantization using that many levels.
     * Negative values mean that we should quantize using the recalibration report's quantization level.
     */
    @Argument(fullName="quantize_quals", shortName = "qq", doc = "Quantize quality scores to a given number of levels (with -BQSR)", required=false)
    public int quantizationLevels = 0;

    /**
     * Static quantized quals are entirely separate from the quantize_qual option which uses dynamic binning.
     * The two types of binning should not be used together.
     */
    @Advanced
    @Argument(fullName="static_quantized_quals", shortName = "SQQ", doc = "Use static quantized quality scores to a given number of levels (with -BQSR)", required=false,  exclusiveOf = "quantize_quals", minValue = QualityUtils.MIN_USABLE_Q_SCORE, maxValue = QualityUtils.MAX_QUAL)
    public List<Integer> staticQuantizationQuals = null;

    /**
     * Round down quantized only works with the static_quantized_quals option, and should not be used with
     * the dynamic binning option provided by quantize_quals.  When roundDown = false, rounding is done in
     * probability space to the nearest bin.  When roundDown = true, the value is rounded to the nearest bin
     * that is smaller than the current bin.
     */
    @Hidden
    @Argument(fullName="round_down_quantized", shortName = "RDQ", doc = "Round quals down to nearest quantized qual", required=false, exclusiveOf="quantize_quals")
    public boolean roundDown = false;

    /**
     * Turns off printing of the base insertion and base deletion tags when using the -BQSR argument. Only the base substitution qualities will be produced.
     */
    @Argument(fullName="disable_indel_quals", shortName = "DIQ", doc = "Disable printing of base insertion and deletion tags (with -BQSR)", required=false)
    public boolean disableIndelQuals = false;

    /**
     * By default, the OQ tag in not emitted when using the -BQSR argument. Use this flag to include OQ tags in the output BAM or CRAM file.
     * Note that this may results in significant file size increase.
     */
    @Argument(fullName="emit_original_quals", shortName = "EOQ", doc = "Emit the OQ tag with the original base qualities (with -BQSR)", required=false)
    public boolean emitOriginalQuals = false;

    /**
     * This flag tells GATK not to modify quality scores less than this value. Instead they will be written out unmodified in the recalibrated BAM or CRAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Advanced
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Don't recalibrate bases with quality scores less than this threshold (with -BQSR)", required = false, minValue = 0, minRecommendedValue = QualityUtils.MIN_USABLE_Q_SCORE)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;
    /**
     * If specified, this value will be used as the prior for all mismatch quality scores instead of the actual reported quality score.
     */
    @Advanced
    @Argument(fullName = "globalQScorePrior", shortName = "globalQScorePrior", doc = "Global Qscore Bayesian prior to use for BQSR", required = false)
    public double globalQScorePrior = -1.0;


    // --------------------------------------------------------------------------------------------------------------
    //
    // Other utility arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Keep in mind that if you set this to LENIENT, we may refuse to provide you with support if anything goes wrong.
     */
    @Argument(fullName = "validation_strictness", shortName = "S", doc = "How strict should we be with validation", required = false)
    public ValidationStringency strictnessLevel = ValidationStringency.SILENT;
    /**
     * Some tools keep program records in the SAM header by default. Use this argument to override that behavior and discard program records for the SAM header. Does not work on CRAM files.
     */
    @Argument(fullName = "remove_program_records", shortName = "rpr", doc = "Remove program records from the SAM header", required = false)
    public boolean removeProgramRecords = false;
    /**
     * Some tools discard program records from the SAM header by default. Use this argument to override that behavior and keep program records in the SAM header.
     */
    @Argument(fullName = "keep_program_records", shortName = "kpr", doc = "Keep program records in the SAM header", required = false)
    public boolean keepProgramRecords = false;

    /**
     * On-the-fly sample renaming works only with single-sample BAM, CRAM, and VCF files. Each line of the mapping file
     * must contain the absolute path to a BAM, CRAM, or VCF file, followed by whitespace, followed by the new sample
     * name for that BAM, CRAM, or VCF file. The sample name may contain non-tab whitespace, but leading or trailing
     * whitespace will be ignored. The engine will verify at runtime that each BAM/CRAM/VCF targeted for sample
     * renaming has only a single sample specified in its header (though, in the case of BAM/CRAM files, there may be
     * multiple read groups for that sample).
     */
    @Advanced
    @Argument(fullName = "sample_rename_mapping_file", shortName = "sample_rename_mapping_file", doc = "Rename sample IDs on-the-fly at runtime using the provided mapping file", required = false)
    public File sampleRenameMappingFile = null;

    /**
     * For expert users only who know what they are doing. We do not support usage of this argument, so we may refuse to help you if you use it and something goes wrong. The one exception to this rule is ALLOW_N_CIGAR_READS, which is necessary for RNAseq analysis.
     */
    @Advanced
    @Argument(fullName = "unsafe", shortName = "U", doc = "Enable unsafe operations: nothing will be checked at runtime", required = false)
    public ValidationExclusion.TYPE unsafe;
    /**
     * Not recommended for general use. Disables both auto-generation of index files and index file locking
     * when reading VCFs and other rods and an index isn't present or is out-of-date. The file locking necessary for auto index
     * generation to work safely is prone to random failures/hangs on certain platforms, which makes it desirable to disable it
     * for situations like test suite runs where the indices are already known to exist, however this option is unsafe in general
     * because it allows reading from index files without first acquiring a lock.
     */
    @Advanced
    @Argument(fullName = "disable_auto_index_creation_and_locking_when_reading_rods", shortName = "disable_auto_index_creation_and_locking_when_reading_rods",
              doc = "Disable both auto-generation of index files and index file locking",
              required = false)
    public boolean disableAutoIndexCreationAndLockingWhenReadingRods = false;

    @Hidden
    @Argument(fullName = "no_cmdline_in_header", shortName = "no_cmdline_in_header", doc = "Don't output the usual VCF header tag with the command line. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.",
              required = false)
    public boolean disableCommandLineInVCF = false;

    @Argument(fullName = "sites_only", shortName = "sites_only", doc = "Just output sites without genotypes (i.e. only the first 8 columns of the VCF)",
              required = false)
    public boolean sitesOnlyVCF = false;

    /**
     * <p>The VCF specification permits missing records to be dropped from the end of FORMAT fields, so long as GT is always output.
     * This option prevents GATK from performing that trimming.</p>
     *
     * <p>For example, given a FORMAT of <pre>GT:AD:DP:PL</pre>, GATK will by default emit <pre>./.</pre> for a variant with
     * no reads present (ie, the AD, DP, and PL fields are trimmed).  If you specify -writeFullFormat, this record
     * would be emitted as <pre>./.:.:.:.</pre></p>
     */
    @Argument(fullName = "never_trim_vcf_format_field", shortName = "writeFullFormat", doc = "Always output all the records in VCF FORMAT fields, even if some are missing",
              required = false)
    public boolean neverTrimVCFFormatField = false;

    @Hidden
    @Argument(fullName = "bcf", shortName = "bcf", doc = "Force BCF output, regardless of the file's extension",
              required = false)
    public boolean forceBCFOutput = false;

    @Advanced
    @Argument(fullName = "bam_compression", shortName = "compress", doc = "Compression level to use for writing BAM files (0 - 9, higher is more compressed)",
              minValue = 0, maxValue = 9, required = false)
    public Integer bamCompression = null;

    @Advanced
    @Argument(fullName = "simplifyBAM", shortName = "simplifyBAM",
              doc = "If provided, output BAM/CRAM files will be simplified to include just key reads for downstream variation discovery analyses (removing duplicates, PF-, non-primary reads), as well stripping all extended tags from the kept reads except the read group identifier",
              required = false)
    public boolean simplifyBAM = false;

    @Advanced
    @Argument(fullName = "disable_bam_indexing", doc = "Turn off on-the-fly creation of indices for output BAM/CRAM files.",
            required = false)
    public boolean disableBAMIndexing = false;

    @Argument(fullName = "generate_md5", doc = "Enable on-the-fly creation of md5s for output BAM files.",
            required = false)
    public boolean enableBAMmd5 = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Multi-threading arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Data threads contains N cpu threads per data thread, and act as completely data parallel processing, increasing
     * the memory usage of GATK by M data threads. Data threads generally scale extremely effectively, up to 24 cores.
     * See online documentation FAQs for more information.
     */
    @Argument(fullName = "num_threads", shortName = "nt", doc = "Number of data threads to allocate to this analysis", required = false, minValue = 1)
    public Integer numberOfDataThreads = 1;

    /**
     * Each CPU thread operates the map cycle independently, but may run into earlier scaling problems with IO than
     * data threads. Has the benefit of not requiring X times as much memory per thread as data threads do, but rather
     * only a constant overhead. See online documentation FAQs for more information.
     */
    @Argument(fullName="num_cpu_threads_per_data_thread", shortName = "nct", doc="Number of CPU threads to allocate per data thread", required = false, minValue = 1)
    public int numberOfCPUThreadsPerDataThread = 1;

    @Argument(fullName="num_io_threads", shortName = "nit", doc="Number of given threads to allocate to BAM IO", required = false, minValue = 0)
    @Hidden
    public int numberOfIOThreads = 0;

    /**
     * Enable GATK to monitor its own threading efficiency, at an itsy-bitsy tiny
     * cost (< 0.1%) in runtime because of turning on the JavaBean.  This is largely for
     * debugging purposes. Note that this argument is not compatible with -nt, it only works with -nct.
     */
    @Argument(fullName = "monitorThreadEfficiency", shortName = "mte", doc = "Enable threading efficiency monitoring", required = false)
    public Boolean monitorThreadEfficiency = false;

    @Argument(fullName = "num_bam_file_handles", shortName = "bfh", doc="When using IO threads, total number of BAM file handles to keep open simultaneously", required=false, minValue = 1)
    @Hidden
    public Integer numberOfBAMFileHandles = null;
    /**
     * This will filter out read groups matching <TAG>:<STRING> (e.g. SM:sample1) or a .txt file containing the filter strings one per line.
     */
    @Input(fullName = "read_group_black_list", shortName="rgbl", doc="Exclude read groups based on tags", required = false)
    public List<String> readGroupBlackList = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // PED (pedigree) support
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * <p>Reads PED file-formatted tabular text files describing meta-data about the samples being
     * processed in the GATK.</p>
     *
     * <ul>
     *  <li>see <a href="http://www.broadinstitute.org/mpg/tagger/faq.html">http://www.broadinstitute.org/mpg/tagger/faq.html</a></li>
     *  <li>see <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped">http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped</a></li>
     * </ul>
     *
     * <p>The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:</p>
     *
     * <ul>
     *  <li>Family ID</li>
     *  <li>Individual ID</li>
     *  <li>Paternal ID</li>
     *  <li>Maternal ID</li>
     *  <li>Sex (1=male; 2=female; other=unknown)</li>
     *  <li>Phenotype</li>
     * </ul>
     *
     *  <p>The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a person.
     *  A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a
     *  quantitative trait or an affection status column: GATK will automatically detect which type
     *  (i.e. based on whether a value other than 0, 1, 2 or the missing genotype code is observed).</p>
     *
     *  <p>If an individual's sex is unknown, then any character other than 1 or 2 can be used.</p>
     *
     *  <p>You can add a comment to a PED or MAP file by starting the line with a # character. The rest of that
     *  line will be ignored. Do not start any family IDs with this character therefore.</p>
     *
     *  <p>Affection status should be coded:</p>
     *
     * <ul>
     *  <li>-9 missing</li>
     *   <li>0 missing</li>
     *   <li>1 unaffected</li>
     *   <li>2 affected</li>
     * </ul>
     *
     * <p>If any value outside of -9,0,1,2 is detected than the samples are assumed
     * to phenotype values are interpreted as string phenotype values.  In this case -9 uniquely
     * represents the missing value.</p>
     *
     * <p>Genotypes (column 7 onwards) cannot be specified to the GATK.</p>
     *
     * <p>For example, here are two individuals (one row = one person):</p>
     *
     * <pre>
     *   FAM001  1  0 0  1  2
     *   FAM001  2  0 0  1  2
     * </pre>
     *
     * <p>Each -ped argument can be tagged with NO_FAMILY_ID, NO_PARENTS, NO_SEX, NO_PHENOTYPE to
     * tell the GATK PED parser that the corresponding fields are missing from the ped file.</p>
     *
     * <p>Note that most GATK walkers do not use pedigree information.  Walkers that require pedigree
     * data should clearly indicate so in their arguments and will throw errors if required pedigree
     * information is missing.</p>
     */
    @Argument(fullName="pedigree", shortName = "ped", doc="Pedigree files for samples",required=false)
    public List<File> pedigreeFiles = Collections.emptyList();

    /**
     * Inline PED records (see -ped argument).  Each -pedString STRING can contain one or more
     * valid PED records (see -ped) separated by semi-colons.  Supports all tags for each pedString
     * as -ped supports
     */
    @Argument(fullName="pedigreeString", shortName = "pedString", doc="Pedigree string for samples",required=false)
    public List<String> pedigreeStrings = Collections.emptyList();

    /**
     * How strict should we be in parsing the PED files?
     */
    @Argument(fullName="pedigreeValidationType", shortName = "pedValidationType", doc="Validation strictness for pedigree information",required=false)
    public PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    // --------------------------------------------------------------------------------------------------------------
    //
    // BAM indexing and sharding arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * NO INTEGRATION TESTS are available.  Use at your own risk.
     */
    @Argument(fullName="allow_intervals_with_unindexed_bam",doc="Allow interval processing with an unsupported BAM/CRAM",required=false)
    @Hidden
    public boolean allowIntervalsWithUnindexedBAM = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing BCF2
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * If provided, whenever we create a VCFWriter we will also write out a BCF file alongside it, for testing purposes.
     */
    @Argument(fullName="generateShadowBCF",shortName = "generateShadowBCF",doc="Write a BCF copy of the output VCF",required=false)
    @Hidden
    public boolean generateShadowBCF = false;
    // TODO -- remove all code tagged with TODO -- remove me when argument generateShadowBCF is removed

    // --------------------------------------------------------------------------------------------------------------
    //
    // VCF/BCF index parameters
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Specify the Tribble indexing strategy to use for VCFs.
     *
     * LINEAR creates a LinearIndex with bins of equal width, specified by the Bin Width parameter
     * INTERVAL creates an IntervalTreeIndex with bins with an equal amount of features, specified by the Features Per Bin parameter
     * DYNAMIC_SEEK attempts to optimize for minimal seek time by choosing an appropriate strategy and parameter (user-supplied parameter is ignored)
     * DYNAMIC_SIZE attempts to optimize for minimal index size by choosing an appropriate strategy and parameter (user-supplied parameter is ignored)
     *
     * This argument is deprecated, using the output file ".g.vcf" extension will automatically set the appropriate value
     */
    @Argument(fullName="variant_index_type",shortName = "variant_index_type",doc="Type of IndexCreator to use for VCF/BCF indices",required=false)
    @Advanced
    public GATKVCFIndexType variant_index_type = GATKVCFUtils.DEFAULT_INDEX_TYPE;
    /**
     * This is either the bin width or the number of features per bin, depending on the indexing strategy
     *
     * This argument is deprecated, using the output file ".g.vcf" extension will automatically set the appropriate value
     */
    @Argument(fullName="variant_index_parameter",shortName = "variant_index_parameter",doc="Parameter to pass to the VCF/BCF IndexCreator",required=false)
    @Advanced
    public int variant_index_parameter = GATKVCFUtils.DEFAULT_INDEX_PARAMETER;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Window arguments
    //
    // -------------------------------------------------------------------------------------------------------------
    /**
     * Stop of the expanded window for which the reference context should be provided, relative to the locus.
     */
    @Argument(fullName = "reference_window_stop", shortName = "ref_win_stop", doc = "Reference window stop", minValue = 0, required = false)
    @Advanced
    public int reference_window_stop = DEFAULT_REFERENCE_WINDOW_STOP;
}

