package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.simpleframework.xml.*;
import org.simpleframework.xml.core.Persister;
import org.simpleframework.xml.stream.Format;
import org.simpleframework.xml.stream.HyphenStyle;

import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * User: aaron
 * Date: May 7, 2009
 * Time: 11:46:21 AM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 */
@Root
public class GATKArgumentCollection {

    /* our version number */
    private float versionNumber = 1;
    private String description = "GATK Arguments";

    /** the constructor */
    public GATKArgumentCollection() {
    }

    @ElementMap(entry = "analysis_argument", key = "key", attribute = true, inline = true, required = false)
    public Map<String, String> walkerArgs = new HashMap<String, String>();

    // parameters and their defaults
    @ElementList(required = false)
    @Argument(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = false)
    public List<File> samFiles = new ArrayList<File>();

    @ElementList(required = false)
    @Argument(fullName = "read_filter", shortName = "rf", doc = "Specify filtration criteria to apply to each read individually.", required = false)
    public List<String> readFilters = new ArrayList<String>();

    @ElementList(required = false)
    @Argument(fullName = "intervals", shortName = "L", doc = "A list of genomic intervals over which to operate. Can be explicitly specified on the command line or in a file.", required = false)
    public List<String> intervals = null;

    @Element(required = false)
    @Argument(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = false)
    public File referenceFile = null;

    @ElementList(required = false)
    @Argument(fullName = "rodBind", shortName = "B", doc = "Bindings for reference-ordered data, in the form <name>,<type>,<file>", required = false)
    public ArrayList<String> RODBindings = new ArrayList<String>();

    @Element(required = false)
    @Argument(fullName = "DBSNP", shortName = "D", doc = "DBSNP file", required = false)
    public String DBSNPFile = null;

    @Element(required = false)
    @Argument(fullName = "hapmap", shortName = "H", doc = "Hapmap file", required = false)
    public String HAPMAPFile = null;

    @Element(required = false)
    @Argument(fullName = "hapmap_chip", shortName = "hc", doc = "Hapmap chip file", required = false)
    public String HAPMAPChipFile = null;

    /** An output file presented to the walker. */
    @Element(required = false)
    @Argument(fullName = "out", shortName = "o", doc = "An output file presented to the walker.  Will overwrite contents if file exists.", required = false)
    public String outFileName = null;

    /** An error output file presented to the walker. */
    @Element(required = false)
    @Argument(fullName = "err", shortName = "e", doc = "An error output file presented to the walker.  Will overwrite contents if file exists.", required = false)
    public String errFileName = null;

    /** A joint file for both 'normal' and error output presented to the walker. */
    @Element(required = false)
    @Argument(fullName = "outerr", shortName = "oe", doc = "A joint file for 'normal' and error output presented to the walker.  Will overwrite contents if file exists.", required = false)
    public String outErrFileName = null;

    @Element(required = false)
    @Argument(fullName = "maximum_iterations", shortName = "M", doc = "Maximum number of iterations to process before exiting, the lower bound is zero.  Intended only for testing", required = false)
    public Integer maximumEngineIterations = -1;

    @Element(required = false)
    @Argument(fullName = "filterZeroMappingQualityReads", shortName = "fmq0", doc = "If true, mapping quality zero reads will be filtered at the lowest GATK level.  Vastly improves performance at areas with abnormal depth due to mapping Q0 reads", required = false)
    public Boolean filterZeroMappingQualityReads = false;

    @Element(required = false)
    @Argument(fullName = "downsample_to_fraction", shortName = "dfrac", doc = "Fraction [0.0-1.0] of reads to downsample to", required = false)
    public Double downsampleFraction = null;

    @Element(required = false)
    @Argument(fullName = "downsample_to_coverage", shortName = "dcov", doc = "Coverage [integer] to downsample to", required = false)
    public Integer downsampleCoverage = null;

    @Element(required = false)
    @Argument(fullName = "validation_strictness", shortName = "S", doc = "How strict should we be with validation (LENIENT|SILENT|STRICT)", required = false)
    public SAMFileReader.ValidationStringency strictnessLevel = SAMFileReader.ValidationStringency.SILENT;

    @Element(required = false)
    @Argument(fullName = "unsafe", shortName = "U", doc = "If set, enables unsafe operations, nothing will be checked at runtime.", required = false)
    public Boolean unsafe = false;

    @Element(required = false)
    @Argument(fullName = "max_reads_at_locus", shortName = "mrl", doc = "Sets the upper limit for the number of reads presented at a single locus. int.MAX_VALUE by default.", required = false)
    public int readMaxPileup = Integer.MAX_VALUE;

    @Element(required = false)
    @Argument(fullName = "disablethreading", shortName = "dt", doc = "Disable experimental threading support.", required = false)
    public Boolean disableThreading = false;

    /** How many threads should be allocated to this analysis. */
    @Element(required = false)
    @Argument(fullName = "numthreads", shortName = "nt", doc = "How many threads should be allocated to running this analysis.", required = false)
    public int numberOfThreads = 1;

    /** What rule should we use when merging intervals */
    @Element(required = false)
    @Argument(fullName = "interval_merging", shortName = "im", doc = "What interval merging rule should we use {ALL [DEFAULT],OVERLAPPING_ONLY,NONE}.", required = false)
    public INTERVAL_MERGING_RULE intervalMerging = INTERVAL_MERGING_RULE.ALL;

    /** Should we enable rodWalkers?  This is currently unsafe */
    @Element(required = false)
    @Argument(fullName = "enableRodWalkers", shortName = "erw", doc = "Enable experimental rodWalker support.  TEMPORARY HACK TO ALLOW EXPERIMENTATION WITH ROD WALKERS.  [default is false]}.", required = false)
    public boolean enableRodWalkers = false;

    /**
     * marshal the data out to a object
     *
     * @param collection the GATKArgumentCollection to load into
     * @param outputFile the file to write to
     */
    public static void marshal(GATKArgumentCollection collection, String outputFile) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        File result = new File(outputFile);
        try {
            serializer.write(collection, result);
        } catch (Exception e) {
            throw new StingException("Failed to marshal the data to the file " + outputFile, e);
        }
    }

    /**
     * marshal the data out to a object
     *
     * @param collection the GATKArgumentCollection to load into
     * @param outputFile the stream to write to
     */
    public static void marshal(GATKArgumentCollection collection, PrintStream outputFile) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        try {
            serializer.write(collection, outputFile);
        } catch (Exception e) {
            throw new StingException("Failed to marshal the data to the file " + outputFile, e);
        }
    }

    /**
     * unmashall the object from a configuration file
     *
     * @param filename the filename to marshal from
     */
    public static GATKArgumentCollection unmarshal(String filename) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        File source = new File(filename);
        try {
            GATKArgumentCollection example = serializer.read(GATKArgumentCollection.class, source);
            return example;
        } catch (Exception e) {
            throw new StingException("Failed to marshal the data from file " + filename, e);
        }
    }

    /**
     * unmashall the object from a configuration file
     *
     * @param file the inputstream to marshal from
     */
    public static GATKArgumentCollection unmarshal(InputStream file) {
        Serializer serializer = new Persister(new Format(new HyphenStyle()));
        try {
            GATKArgumentCollection example = serializer.read(GATKArgumentCollection.class, file);
            return example;
        } catch (Exception e) {
            throw new StingException("Failed to marshal the data from file " + file.toString(), e);
        }
    }


    /**
     * test equality between two arg collections.  This function defines the statement:
     * "not fun to write"
     *
     * @param other the other collection
     *
     * @return true if they're equal
     */
    public boolean equals(GATKArgumentCollection other) {
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
        if (other.RODBindings.size() != RODBindings.size()) {
            return false;
        }
        for (int x = 0; x < RODBindings.size(); x++) {
            if (!RODBindings.get(x).equals(other.RODBindings.get(x))) {
                return false;
            }
        }
        if (!other.samFiles.equals(this.samFiles)) {
            return false;
        }
        if (!other.maximumEngineIterations.equals(this.maximumEngineIterations)) {
            return false;
        }
        if (!other.strictnessLevel.equals(this.strictnessLevel)) {
            return false;
        }
        if (!other.referenceFile.equals(this.referenceFile)) {
            return false;
        }
        if (!other.intervals.equals(this.intervals)) {
            return false;
        }
        if (!other.DBSNPFile.equals(this.DBSNPFile)) {
            return false;
        }
        if (!other.HAPMAPFile.equals(this.HAPMAPFile)) {
            return false;
        }
        if (!other.HAPMAPChipFile.equals(this.HAPMAPChipFile)) {
            return false;
        }
        if (!other.unsafe.equals(this.unsafe)) {
            return false;
        }
        if (other.readMaxPileup != this.readMaxPileup) {
            return false;
        }
        if ((other.filterZeroMappingQualityReads == null && this.filterZeroMappingQualityReads != null) ||
                (other.filterZeroMappingQualityReads != null && !other.filterZeroMappingQualityReads.equals(this.filterZeroMappingQualityReads))) {
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
        if (!other.outFileName.equals(this.outFileName)) {
            return false;
        }
        if (!other.errFileName.equals(this.errFileName)) {
            return false;
        }
        if (!other.outErrFileName.equals(this.outErrFileName)) {
            return false;
        }
        if (other.numberOfThreads != this.numberOfThreads) {
            return false;
        }
        if (other.intervalMerging != this.intervalMerging) {
            return false;
        }
        if (other.enableRodWalkers != this.enableRodWalkers) {
            return false;
        }

        return true;
    }

    /**
     * a class we use to determine the merging rules for intervals passed to the GATK
     */
    public enum INTERVAL_MERGING_RULE {
        ALL, // we merge both overlapping intervals and abutting intervals
        OVERLAPPING_ONLY, // We merge intervals that are overlapping, but NOT ones that only abut each other
        NONE; // we merge neither overlapping or abutting intervals, the list of intervals is sorted, but not merged

        public boolean check() {
            if (this.compareTo(NONE) == 0)
                throw new UnsupportedOperationException("We Currently do not support INTERVAL_MERGING_RULE.NONE");
            return true;
        }
    }
}

