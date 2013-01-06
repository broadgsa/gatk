/*
 * Copyright (c) 2010 The Broad Institute
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
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.gatk.samples.PedigreeValidationType;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;

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

    @Argument(fullName = "phone_home", shortName = "et", doc="What kind of GATK run report should we generate? STANDARD is the default, can be NO_ET so nothing is posted to the run repository. Please see " + GATKRunReport.PHONE_HOME_DOCS_URL + " for details.", required = false)
    public GATKRunReport.PhoneHomeOption phoneHomeType = GATKRunReport.PhoneHomeOption.STANDARD;

    @Argument(fullName = "gatk_key", shortName = "K", doc="GATK Key file. Required if running with -et NO_ET. Please see " + GATKRunReport.PHONE_HOME_DOCS_URL + " for details.", required = false)
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

    /**
     * Using this option one can instruct the GATK engine to traverse over only part of the genome.  This argument can be specified multiple times.
     * One may use samtools-style intervals either explicitly (e.g. -L chr1 or -L chr1:100-200) or listed in a file (e.g. -L myFile.intervals).
     * Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).
     * To specify the completely unmapped reads in the BAM file (i.e. those without a reference contig) use -L unmapped.
     */
    @Input(fullName = "intervals", shortName = "L", doc = "One or more genomic intervals over which to operate. Can be explicitly specified on the command line or in a file (including a rod file)", required = false)
    public List<IntervalBinding<Feature>> intervals = null;

    /**
     * Using this option one can instruct the GATK engine NOT to traverse over certain parts of the genome.  This argument can be specified multiple times.
     * One may use samtools-style intervals either explicitly (e.g. -XL chr1 or -XL chr1:100-200) or listed in a file (e.g. -XL myFile.intervals).
     * Additionally, one may specify a rod file to skip over the positions for which there is a record in the file (e.g. -XL file.vcf).
     */
    @Input(fullName = "excludeIntervals", shortName = "XL", doc = "One or more genomic intervals to exclude from processing. Can be explicitly specified on the command line or in a file (including a rod file)", required = false)
    public List<IntervalBinding<Feature>> excludeIntervals = null;

    /**
     * How should the intervals specified by multiple -L or -XL arguments be combined?  Using this argument one can, for example, traverse over all of the positions
     * for which there is a record in a VCF but just in chromosome 20 (-L chr20 -L file.vcf -isr INTERSECTION).
     */
    @Argument(fullName = "interval_set_rule", shortName = "isr", doc = "Indicates the set merging approach the interval parser should use to combine the various -L or -XL inputs", required = false)
    public IntervalSetRule intervalSetRule = IntervalSetRule.UNION;

    /**
     * Should abutting (but not overlapping) intervals be treated as separate intervals?
     */
    @Argument(fullName = "interval_merging", shortName = "im", doc = "Indicates the interval merging rule we should use for abutting intervals", required = false)
    public IntervalMergingRule intervalMerging = IntervalMergingRule.ALL;

    /**
     * For example, '-L chr1:100' with a padding value of 20 would turn into '-L chr1:80-120'.
     */
    @Argument(fullName = "interval_padding", shortName = "ip", doc = "Indicates how many basepairs of padding to include around each of the intervals specified with the -L/--intervals argument", required = false)
    public int intervalPadding = 0;

    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;

    @Argument(fullName = "nonDeterministicRandomSeed", shortName = "ndrs", doc = "Makes the GATK behave non deterministically, that is, the random numbers generated will be different in every run", required = false)
    public boolean nonDeterministicRandomSeed = false;

    @Argument(fullName = "disableRandomization",doc="Completely eliminates randomization from nondeterministic methods. To be used mostly in the testing framework where dynamic parallelism can result in differing numbers of calls to the generator.")
    public boolean disableRandomization = false;

    @Argument(fullName = "maxRuntime", shortName = "maxRuntime", doc="If provided, that GATK will stop execution cleanly as soon after maxRuntime has been exceeded, truncating the run but not exiting with a failure.  By default the value is interpreted in minutes, but this can be changed by maxRuntimeUnits", required = false)
    public long maxRuntime = GenomeAnalysisEngine.NO_RUNTIME_LIMIT;

    @Argument(fullName = "maxRuntimeUnits", shortName = "maxRuntimeUnits", doc="The TimeUnit for maxRuntime", required = false)
    public TimeUnit maxRuntimeUnits = TimeUnit.MINUTES;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Downsampling Arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    @Argument(fullName = "downsampling_type", shortName="dt", doc="Type of reads downsampling to employ at a given locus.  Reads will be selected randomly to be removed from the pile based on the method described here", required = false)
    public DownsampleType downsamplingType = null;

    @Argument(fullName = "downsample_to_fraction", shortName = "dfrac", doc = "Fraction [0.0-1.0] of reads to downsample to", required = false)
    public Double downsampleFraction = null;

    @Argument(fullName = "downsample_to_coverage", shortName = "dcov", doc = "Coverage [integer] to downsample to at any given locus; note that downsampled reads are randomly selected from all possible reads at a locus. For non-locus-based traversals (eg., ReadWalkers), this sets the maximum number of reads at each alignment start position.", required = false)
    public Integer downsampleCoverage = null;

    @Argument(fullName = "use_legacy_downsampler", shortName = "use_legacy_downsampler", doc = "Use the legacy downsampling implementation instead of the newer, less-tested implementation", required = false)
    public boolean useLegacyDownsampler = false;

    /**
     * Gets the downsampling method explicitly specified by the user.  If the user didn't specify
     * a default downsampling mechanism, return the default.
     * @return The explicitly specified downsampling mechanism, or the default if none exists.
     */
    public DownsamplingMethod getDownsamplingMethod() {
        if ( downsamplingType == null && downsampleFraction == null && downsampleCoverage == null )
            return null;

        return new DownsamplingMethod(downsamplingType, downsampleCoverage, downsampleFraction, useLegacyDownsampler);
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
        useLegacyDownsampler = method.useLegacyDownsampler;
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

    @Argument(fullName = "allow_potentially_misencoded_quality_scores", shortName="allowPotentiallyMisencodedQuals", doc="Do not fail when encountered base qualities that are too high and seemingly indicate a problem with the base quality encoding of the BAM file", required = false)
    public boolean ALLOW_POTENTIALLY_MISENCODED_QUALS = false;

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

    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "If set, use the original base quality scores from the OQ tag when present instead of the standard scores", required=false)
    public Boolean useOriginalBaseQualities = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // BQSR arguments
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Enables on-the-fly recalibrate of base qualities.  The covariates tables are produced by the BaseQualityScoreRecalibrator tool.
     * Please be aware that one should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Input(fullName="BQSR", shortName="BQSR", required=false, doc="The input covariates table file which enables on-the-fly base quality score recalibration")
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

    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "If reads are missing some or all base quality scores, this value will be used for all base quality scores", required=false)
    public byte defaultBaseQualities = -1;

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

    @Argument(fullName = "unsafe", shortName = "U", doc = "If set, enables unsafe operations: nothing will be checked at runtime.  For expert users only who know what they are doing.  We do not support usage of this argument.", required = false)
    public ValidationExclusion.TYPE unsafe;

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
     * Enable GATK to monitor its own threading efficiency, at a itsy-bitsy tiny
     * cost (< 0.1%) in runtime because of turning on the JavaBean.  This is largely for
     * debugging purposes.
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
}

