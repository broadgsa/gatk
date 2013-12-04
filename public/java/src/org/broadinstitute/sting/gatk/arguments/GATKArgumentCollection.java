/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.arguments;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.gatk.samples.PedigreeValidationType;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variant.GATKVCFIndexType;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;

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

    /* our version number */
    private float versionNumber = 1;
    private String description = "GATK Arguments";

    /** the constructor */
    public GATKArgumentCollection() {
    }

    // parameters and their defaults
    @Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = false)
    public List<String> samFiles = new ArrayList<String>();

    @Argument(fullName = "read_buffer_size", shortName = "rbs", doc="Number of reads per SAM file to buffer in memory", required = false)
    public Integer readBufferSize = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // GATKRunReport options
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(fullName = "phone_home", shortName = "et", doc="What kind of GATK run report should we generate? AWS is the default, can be NO_ET so nothing is posted to the run repository. Please see " + UserException.PHONE_HOME_DOCS_URL + " for details.", required = false)
    public GATKRunReport.PhoneHomeOption phoneHomeType = GATKRunReport.PhoneHomeOption.AWS;

    @Argument(fullName = "gatk_key", shortName = "K", doc="GATK Key file. Required if running with -et NO_ET. Please see " + UserException.PHONE_HOME_DOCS_URL + " for details.", required = false)
    public File gatkKeyFile = null;

    /**
     * The GATKRunReport supports (as of GATK 2.2) tagging GATK runs with an arbitrary String tag that can be
     * used to group together runs during later analysis.  One use of this capability is to tag runs as GATK
     * performance tests, so that the performance of the GATK over time can be assessed from the logs directly.
     *
     * Note that the tags do not conform to any ontology, so you are free to use any tags that you might find
     * meaningful.
     */
    @Argument(fullName = "tag", shortName = "tag", doc="Arbitrary tag string to identify this GATK run as part of a group of runs, for later analysis", required = false)
    public String tag = "NA";

    // --------------------------------------------------------------------------------------------------------------
    //
    // General features
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(fullName = "read_filter", shortName = "rf", doc = "Specify filtration criteria to apply to each read individually", required = false)
    public List<String> readFilters = new ArrayList<String>();

    @ArgumentCollection
    public IntervalArgumentCollection intervalArguments = new IntervalArgumentCollection();

    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;

    @Argument(fullName = "nonDeterministicRandomSeed", shortName = "ndrs", doc = "Makes the GATK behave non deterministically, that is, the random numbers generated will be different in every run", required = false)
    public boolean nonDeterministicRandomSeed = false;

    @Hidden
    @Argument(fullName = "disableDithering",doc="Completely eliminates randomized dithering from rank sum tests. To be used in the testing framework where dynamic parallelism can result in differing numbers of calls to the random generator.")
    public boolean disableDithering = false;

    @Argument(fullName = "maxRuntime", shortName = "maxRuntime", doc="If provided, that GATK will stop execution cleanly as soon after maxRuntime has been exceeded, truncating the run but not exiting with a failure.  By default the value is interpreted in minutes, but this can be changed by maxRuntimeUnits", required = false)
    public long maxRuntime = GenomeAnalysisEngine.NO_RUNTIME_LIMIT;

    @Argument(fullName = "maxRuntimeUnits", shortName = "maxRuntimeUnits", doc="The TimeUnit for maxRuntime", required = false)
    public TimeUnit maxRuntimeUnits = TimeUnit.MINUTES;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Downsampling Arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Reads will be selected randomly to be removed from the pile based on the method described here.
     */
    @Argument(fullName = "downsampling_type", shortName="dt", doc="Type of reads downsampling to employ at a given locus", required = false)
    public DownsampleType downsamplingType = null;

    @Argument(fullName = "downsample_to_fraction", shortName = "dfrac", doc = "Fraction [0.0-1.0] of reads to downsample to", required = false)
    public Double downsampleFraction = null;

    /**
     * For locus-based traversals (LocusWalkers and ActiveRegionWalkers), downsample_to_coverage controls the
     * maximum depth of coverage at each locus. For read-based traversals (ReadWalkers), it controls the
     * maximum number of reads sharing the same alignment start position. For ReadWalkers you will typically need to use
     * much lower dcov values than you would with LocusWalkers to see an effect. Note that this downsampling option does
     * not produce an unbiased random sampling from all available reads at each locus: instead, the primary goal of the
     * to-coverage downsampler is to maintain an even representation of reads from all alignment start positions when
     * removing excess coverage. For a truly unbiased random sampling of reads, use -dfrac instead. Also note
     * that the coverage target is an approximate goal that is not guaranteed to be met exactly: the downsampling
     * algorithm will under some circumstances retain slightly more or less coverage than requested.
     */
    @Argument(fullName = "downsample_to_coverage", shortName = "dcov",
              doc = "Coverage [integer] to downsample to per locus (for locus walkers) or per alignment start position (for read walkers)",
              required = false)
    public Integer downsampleCoverage = null;

    /**
     * Gets the downsampling method explicitly specified by the user.  If the user didn't specify
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
    @Argument(fullName = "baq", shortName="baq", doc="Type of BAQ calculation to apply in the engine", required = false)
    public BAQ.CalculationMode BAQMode = BAQ.CalculationMode.OFF;

    @Argument(fullName = "baqGapOpenPenalty", shortName="baqGOP", doc="BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets", required = false)
    public double BAQGOP = BAQ.DEFAULT_GOP;

    // --------------------------------------------------------------------------------------------------------------
    //
    // quality encoding checking arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Q0 == ASCII 33 according to the SAM specification, whereas Illumina encoding starts at Q64.  The idea here is
     * simple: we just iterate over all reads and subtract 31 from every quality score.
     */
    @Argument(fullName = "fix_misencoded_quality_scores", shortName="fixMisencodedQuals", doc="Fix mis-encoded base quality scores", required = false)
    public boolean FIX_MISENCODED_QUALS = false;

    @Argument(fullName = "allow_potentially_misencoded_quality_scores", shortName="allowPotentiallyMisencodedQuals", doc="Do not fail when encountering base qualities that are too high and that seemingly indicate a problem with the base quality encoding of the BAM file", required = false)
    public boolean ALLOW_POTENTIALLY_MISENCODED_QUALS = false;

    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "If set, use the original base quality scores from the OQ tag when present instead of the standard scores", required=false)
    public Boolean useOriginalBaseQualities = false;

    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "If reads are missing some or all base quality scores, this value will be used for all base quality scores", required=false)
    public byte defaultBaseQualities = -1;

    // --------------------------------------------------------------------------------------------------------------
    //
    // performance log arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * The file name for the GATK performance log output, or null if you don't want to generate the
     * detailed performance logging table.  This table is suitable for importing into R or any
     * other analysis software that can read tsv files
     */
    @Argument(fullName = "performanceLog", shortName="PF", doc="If provided, a GATK runtime performance log will be written to this file", required = false)
    public File performanceLog = null;

    // --------------------------------------------------------------------------------------------------------------
    //
    // BQSR arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Enables on-the-fly recalibrate of base qualities.  The covariates tables are produced by the BaseQualityScoreRecalibrator tool.
     * Please be aware that one should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Input(fullName="BQSR", shortName="BQSR", required=false, doc="The input covariates table file which enables on-the-fly base quality score recalibration (intended for use with BaseRecalibrator and PrintReads)")
    public File BQSR_RECAL_FILE = null;

    /**
     * Turns on the base quantization module. It requires a recalibration report (-BQSR).
     *
     * A value of 0 here means "do not quantize".
     * Any value greater than zero will be used to recalculate the quantization using that many levels.
     * Negative values mean that we should quantize using the recalibration report's quantization level.
     */
    @Hidden
    @Argument(fullName="quantize_quals", shortName = "qq", doc = "Quantize quality scores to a given number of levels (with -BQSR)", required=false)
    public int quantizationLevels = 0;

    /**
     * Turns off printing of the base insertion and base deletion tags when using the -BQSR argument and only the base substitution qualities will be produced.
     */
    @Argument(fullName="disable_indel_quals", shortName = "DIQ", doc = "If true, disables printing of base insertion and base deletion tags (with -BQSR)", required=false)
    public boolean disableIndelQuals = false;

    /**
     * By default, the OQ tag in not emitted when using the -BQSR argument.
     */
    @Argument(fullName="emit_original_quals", shortName = "EOQ", doc = "If true, enables printing of the OQ tag with the original base qualities (with -BQSR)", required=false)
    public boolean emitOriginalQuals = false;

    /**
     * Do not modify quality scores less than this value but rather just write them out unmodified in the recalibrated BAM file.
     * In general it's unsafe to change qualities scores below < 6, since base callers use these values to indicate random or bad bases.
     * For example, Illumina writes Q2 bases when the machine has really gone wrong. This would be fine in and of itself,
     * but when you select a subset of these reads based on their ability to align to the reference and their dinucleotide effect,
     * your Q2 bin can be elevated to Q8 or Q10, leading to issues downstream.
     */
    @Argument(fullName = "preserve_qscores_less_than", shortName = "preserveQ", doc = "Bases with quality scores less than this threshold won't be recalibrated (with -BQSR)", required = false)
    public int PRESERVE_QSCORES_LESS_THAN = QualityUtils.MIN_USABLE_Q_SCORE;

    @Argument(fullName = "globalQScorePrior", shortName = "globalQScorePrior", doc = "The global Qscore Bayesian prior to use in the BQSR. If specified, this value will be used as the prior for all mismatch quality scores instead of the actual reported quality score", required = false)
    public double globalQScorePrior = -1.0;

    /**
     * For the sake of your data, please only use this option if you know what you are doing.  It is absolutely not recommended practice
     * to run base quality score recalibration on reduced BAM files.
     */
    @Advanced
    @Argument(fullName = "allow_bqsr_on_reduced_bams_despite_repeated_warnings", shortName="allowBqsrOnReducedBams", doc="Do not fail when running base quality score recalibration on a reduced BAM file even though we highly recommend against it", required = false)
    public boolean ALLOW_BQSR_ON_REDUCED_BAMS = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Other utility arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(fullName = "validation_strictness", shortName = "S", doc = "How strict should we be with validation", required = false)
    public SAMFileReader.ValidationStringency strictnessLevel = SAMFileReader.ValidationStringency.SILENT;

    @Argument(fullName = "remove_program_records", shortName = "rpr", doc = "Should we override the Walker's default and remove program records from the SAM header", required = false)
    public boolean removeProgramRecords = false;

    @Argument(fullName = "keep_program_records", shortName = "kpr", doc = "Should we override the Walker's default and keep program records from the SAM header", required = false)
    public boolean keepProgramRecords = false;

    @Advanced
    @Argument(fullName = "sample_rename_mapping_file", shortName = "sample_rename_mapping_file",
              doc = "Rename sample IDs on-the-fly at runtime using the provided mapping file. This option requires that " +
                    "each BAM file listed in the mapping file have only a single sample specified in its header (though there " +
                    "may be multiple read groups for that sample). Each line of the mapping file must contain the absolute path " +
                    "to a BAM file, followed by whitespace, followed by the new sample name for that BAM file.",
              required = false)
    public File sampleRenameMappingFile = null;

    @Argument(fullName = "unsafe", shortName = "U", doc = "If set, enables unsafe operations: nothing will be checked at runtime.  For expert users only who know what they are doing.  We do not support usage of this argument.", required = false)
    public ValidationExclusion.TYPE unsafe;

    @Hidden
    @Advanced
    @Argument(fullName = "disable_auto_index_creation_and_locking_when_reading_rods", shortName = "disable_auto_index_creation_and_locking_when_reading_rods",
              doc = "UNSAFE FOR GENERAL USE (FOR TEST SUITE USE ONLY). Disable both auto-generation of index files and index file locking " +
                    "when reading VCFs and other rods and an index isn't present or is out-of-date. The file locking necessary for auto index " +
                    "generation to work safely is prone to random failures/hangs on certain platforms, which makes it desirable to disable it " +
                    "for situations like test suite runs where the indices are already known to exist, however this option is unsafe in general " +
                    "because it allows reading from index files without first acquiring a lock.",
              required = false)
    public boolean disableAutoIndexCreationAndLockingWhenReadingRods = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Multi-threading arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * How many data threads should be allocated to this analysis?  Data threads contains N cpu threads per
     * data thread, and act as completely data parallel processing, increasing the memory usage of GATK
     * by M data threads.  Data threads generally scale extremely effectively, up to 24 cores
     */
    @Argument(fullName = "num_threads", shortName = "nt", doc = "How many data threads should be allocated to running this analysis.", required = false)
    public Integer numberOfDataThreads = 1;

    /**
     * How many CPU threads should be allocated per data thread?  Each CPU thread operates the map
     * cycle independently, but may run into earlier scaling problems with IO than data threads.  Has
     * the benefit of not requiring X times as much memory per thread as data threads do, but rather
     * only a constant overhead.
     */
    @Argument(fullName="num_cpu_threads_per_data_thread", shortName = "nct", doc="How many CPU threads should be allocated per data thread to running this analysis?", required = false)
    public int numberOfCPUThreadsPerDataThread = 1;

    @Argument(fullName="num_io_threads", shortName = "nit", doc="How many of the given threads should be allocated to IO", required = false)
    @Hidden
    public int numberOfIOThreads = 0;

    /**
     * Enable GATK to monitor its own threading efficiency, at an itsy-bitsy tiny
     * cost (< 0.1%) in runtime because of turning on the JavaBean.  This is largely for
     * debugging purposes. Note that this argument is not compatible with -nt, it only works with -nct.
     */
    @Argument(fullName = "monitorThreadEfficiency", shortName = "mte", doc = "Enable GATK threading efficiency monitoring", required = false)
    public Boolean monitorThreadEfficiency = false;

    @Argument(fullName = "num_bam_file_handles", shortName = "bfh", doc="The total number of BAM file handles to keep open simultaneously", required=false)
    public Integer numberOfBAMFileHandles = null;

    @Input(fullName = "read_group_black_list", shortName="rgbl", doc="Filters out read groups matching <TAG>:<STRING> or a .txt file containing the filter strings one per line.", required = false)
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
    @Argument(fullName="pedigreeValidationType", shortName = "pedValidationType", doc="How strict should we be in validating the pedigree information?",required=false)
    public PedigreeValidationType pedigreeValidationType = PedigreeValidationType.STRICT;

    // --------------------------------------------------------------------------------------------------------------
    //
    // BAM indexing and sharding arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(fullName="allow_intervals_with_unindexed_bam",doc="Allow interval processing with an unsupported BAM.  NO INTEGRATION TESTS are available.  Use at your own risk.",required=false)
    @Hidden
    public boolean allowIntervalsWithUnindexedBAM = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing BCF2
    //
    // --------------------------------------------------------------------------------------------------------------

    @Argument(fullName="generateShadowBCF",shortName = "generateShadowBCF",doc="If provided, whenever we create a VCFWriter we will also write out a BCF file alongside it, for testing purposes",required=false)
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
     */

    @Argument(fullName="variant_index_type",shortName = "variant_index_type",doc="which type of IndexCreator to use for VCF/BCF indices",required=false)
    @Advanced
    public GATKVCFIndexType variant_index_type = GATKVCFUtils.DEFAULT_INDEX_TYPE;

    @Argument(fullName="variant_index_parameter",shortName = "variant_index_parameter",doc="the parameter (bin width or features per bin) to pass to the VCF/BCF IndexCreator",required=false)
    @Advanced
    public int variant_index_parameter = GATKVCFUtils.DEFAULT_INDEX_PARAMETER;
}

