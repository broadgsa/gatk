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
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport;
import org.broadinstitute.sting.gatk.samples.PedigreeValidationType;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;

import java.io.File;
import java.util.*;

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

    public Map<String, String> walkerArgs = new HashMap<String, String>();

    // parameters and their defaults
    @Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = false)
    public List<String> samFiles = new ArrayList<String>();

    @Argument(fullName = "read_buffer_size", shortName = "rbs", doc="Number of reads per SAM file to buffer in memory", required = false)
    public Integer readBufferSize = null;

    @Argument(fullName = "phone_home", shortName = "et", doc="What kind of GATK run report should we generate?  Standard is the default, can be verbose or NO_ET so nothing is posted to the run repository", required = false)
    public GATKRunReport.PhoneHomeOption phoneHomeType = GATKRunReport.PhoneHomeOption.STANDARD;

    @Argument(fullName = "read_filter", shortName = "rf", doc = "Specify filtration criteria to apply to each read individually", required = false)
    public List<String> readFilters = new ArrayList<String>();

    /**
     * Using this option one can instruct the GATK engine to traverse over only part of the genome.  This argument can be specified multiple times.
     * One may use samtools-style intervals either explicitly (e.g. -L chr1 or -L chr1:100-200) or listed in a file (e.g. -L myFile.intervals).
     * Additionally, one may specify a rod file to traverse over the positions for which there is a record in the file (e.g. -L file.vcf).
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

    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;

    @Deprecated
    @Hidden
    @Input(fullName = "rodBind", shortName = "B", doc = "Bindings for reference-ordered data, in the form :<name>,<type> <file>", required = false)
    public ArrayList<String> RODBindings = new ArrayList<String>();

    @Argument(fullName = "nonDeterministicRandomSeed", shortName = "ndrs", doc = "Makes the GATK behave non deterministically, that is, the random numbers generated will be different in every run", required = false)
    public boolean nonDeterministicRandomSeed = false;

    /**
     * The override mechanism in the GATK, by default, populates the command-line arguments, then
     * the defaults from the walker annotations.  Unfortunately, walker annotations should be trumped
     * by a user explicitly specifying command-line arguments.
     * TODO: Change the GATK so that walker defaults are loaded first, then command-line arguments.
     */
    private static DownsampleType DEFAULT_DOWNSAMPLING_TYPE = DownsampleType.BY_SAMPLE;
    private static int DEFAULT_DOWNSAMPLING_COVERAGE = 1000;

    @Argument(fullName = "downsampling_type", shortName="dt", doc="Type of reads downsampling to employ at a given locus.  Reads will be selected randomly to be removed from the pile based on the method described here", required = false)
    public DownsampleType downsamplingType = null;

    @Argument(fullName = "downsample_to_fraction", shortName = "dfrac", doc = "Fraction [0.0-1.0] of reads to downsample to", required = false)
    public Double downsampleFraction = null;

    @Argument(fullName = "downsample_to_coverage", shortName = "dcov", doc = "Coverage [integer] to downsample to at any given locus; note that downsampled reads are randomly selected from all possible reads at a locus", required = false)
    public Integer downsampleCoverage = null;

    /**
     * Gets the downsampling method explicitly specified by the user.  If the user didn't specify
     * a default downsampling mechanism, return the default.
     * @return The explicitly specified downsampling mechanism, or the default if none exists.
     */
    public DownsamplingMethod getDownsamplingMethod() {
        if(downsamplingType == null && downsampleFraction == null && downsampleCoverage == null)
            return null;
        if(downsamplingType == null && downsampleCoverage != null)
            return new DownsamplingMethod(DEFAULT_DOWNSAMPLING_TYPE,downsampleCoverage,null);
        return new DownsamplingMethod(downsamplingType,downsampleCoverage,downsampleFraction);
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
    // performance log arguments
    //
    // --------------------------------------------------------------------------------------------------------------
    @Argument(fullName = "performanceLog", shortName="PF", doc="If provided, a GATK runtime performance log will be written to this file", required = false)
    public File performanceLog = null;

    /**
     * Gets the default downsampling method, returned if the user didn't specify any downsampling
     * method.
     * @return The default downsampling mechanism, or null if none exists.
     */
    public static DownsamplingMethod getDefaultDownsamplingMethod() {
        return new DownsamplingMethod(DEFAULT_DOWNSAMPLING_TYPE,DEFAULT_DOWNSAMPLING_COVERAGE,null);
    }

    @Argument(fullName="useOriginalQualities", shortName = "OQ", doc = "If set, use the original base quality scores from the OQ tag when present instead of the standard scores", required=false)
    public Boolean useOriginalBaseQualities = false;

    @Argument(fullName="defaultBaseQualities", shortName = "DBQ", doc = "If reads are missing some or all base quality scores, this value will be used for all base quality scores", required=false)
    public byte defaultBaseQualities = -1;

    @Argument(fullName = "validation_strictness", shortName = "S", doc = "How strict should we be with validation", required = false)
    public SAMFileReader.ValidationStringency strictnessLevel = SAMFileReader.ValidationStringency.SILENT;

    @Argument(fullName = "unsafe", shortName = "U", doc = "If set, enables unsafe operations: nothing will be checked at runtime.  For expert users only who know what they are doing.  We do not support usage of this argument.", required = false)
    public ValidationExclusion.TYPE unsafe;

    @Argument(fullName = "num_threads", shortName = "nt", doc = "How many threads should be allocated to running this analysis", required = false)
    public int numberOfThreads = 1;

    @Input(fullName = "read_group_black_list", shortName="rgbl", doc="Filters out read groups matching <TAG>:<STRING> or a .txt file containing the filter strings one per line", required = false)
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

    @Argument(fullName="disable_experimental_low_memory_sharding",doc="Disable experimental low-memory sharding functionality",required=false)
    public boolean disableLowMemorySharding = false;

    // --------------------------------------------------------------------------------------------------------------
    //
    // methods
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * test equality between two arg collections.  This function defines the statement:
     * "not fun to write"
     *
     * @param other the other collection
     *
     * @return true if they're equal
     */
    public boolean equals(GATKArgumentCollection other) {
        if (other == null) return false;
        if (other.samFiles.size() != samFiles.size()) {
            return false;
        }
        for (int x = 0; x < samFiles.size(); x++) {
            if (!samFiles.get(x).equals(other.samFiles.get(x))) {
                return false;
            }
        }
        if (other.walkerArgs.size() != walkerArgs.size()) {
            return false;
        }
        for (String s : walkerArgs.keySet()) {
            if (!other.walkerArgs.containsKey(s)) {
                return false;
            }
        }
        if (!other.samFiles.equals(this.samFiles)) {
            return false;
        }
        if(other.readBufferSize == null || this.readBufferSize == null) {
            // If either is null, return false if they're both null, otherwise keep going...
            if(other.readBufferSize != null || this.readBufferSize != null)
                return false;
        }
        else {
            if(!other.readBufferSize.equals(this.readBufferSize))
                return false;
        }
        if (!(other.readBufferSize == null && this.readBufferSize == null) && (other.readBufferSize == null || this.readBufferSize == null)) {
            return false;
        }
        if (!other.strictnessLevel.equals(this.strictnessLevel)) {
            return false;
        }
        if (!other.referenceFile.equals(this.referenceFile)) {
            return false;
        }
        if ((other.intervals == null && this.intervals != null) || !other.intervals.equals(this.intervals)) {
            return false;
        }
        if (!other.excludeIntervals.equals(this.excludeIntervals)) {
            return false;
        }
        if (!other.unsafe.equals(this.unsafe)) {
            return false;
        }
        if ((other.downsampleFraction == null && this.downsampleFraction != null) ||
                (other.downsampleFraction != null && !other.downsampleFraction.equals(this.downsampleFraction))) {
            return false;
        }
        if ((other.downsampleCoverage == null && this.downsampleCoverage != null) ||
                (other.downsampleCoverage != null && !other.downsampleCoverage.equals(this.downsampleCoverage))) {
            return false;
        }
        if (other.numberOfThreads != this.numberOfThreads) {
            return false;
        }
        if (other.intervalMerging != this.intervalMerging) {
            return false;
        }

        if (other.phoneHomeType != this.phoneHomeType) {
            return false;
        }

        if (intervalSetRule != other.intervalSetRule)
            return false;

        if ( BAQMode != other.BAQMode ) return false;
        if ( BAQGOP != other.BAQGOP ) return false;

        if ((other.performanceLog == null && this.performanceLog != null) ||
                (other.performanceLog != null && !other.performanceLog.equals(this.performanceLog)))
            return false;

        if (allowIntervalsWithUnindexedBAM != other.allowIntervalsWithUnindexedBAM)
            return false;

        if (disableLowMemorySharding != other.disableLowMemorySharding)
            return false;

        return true;
    }

}

