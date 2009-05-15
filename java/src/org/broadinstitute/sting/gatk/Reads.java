package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
/**
 * User: hanna
 * Date: May 14, 2009
 * Time: 4:06:26 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A data structure containing information about the reads data sources as well as
 * information about how they should be downsampled, sorted, and filtered.
 */
public class Reads {
    private List<File> readsFiles = null;

    private Double downsamplingFraction = null;
    private Integer downsampleToCoverage = null;
    private Integer maxOnFlySorts = null;
    private Boolean beSafe = null;

    /**
     * Gets a list of the files acting as sources of reads.
     * @return A list of files storing reads data.
     */
    public List<File> getReadsFiles() {
        return readsFiles;
    }

    /**
     * Get the fraction of reads to downsample.
     * @return Downsample fraction.
     */
    public Double getDownsamplingFraction() {
        return downsamplingFraction;
    }

    /**
     * Downsample each locus to the specified coverage.
     * @return Coverage to which to downsample.
     */
    public Integer getDownsampleToCoverage() {
        return downsampleToCoverage;
    }

    /**
     * Get the maximum number of supported on-the-fly sorts.
     * @return Max number of on-the-fly sorts.
     */
    public Integer getMaxOnTheFlySorts() {
        return maxOnFlySorts;
    }

    /**
     * Return whether to 'verify' the reads as we pass through them.
     * @return Whether to verify the reads.
     */
    public Boolean getSafetyChecking() {
        return beSafe;
    }

    /**
     * Simple constructor for unit testing.
     * @param readsFiles List of reads files to open.
     */
    public Reads( List<File> readsFiles ) {
        this.readsFiles = readsFiles;
    }

    /**
     * Extract the command-line arguments having to do with reads input
     * files and store them in an easy-to-work-with package.  Constructor
     * is package protected.
     * @param arguments GATK parsed command-line arguments.
     */
    Reads( GATKArgumentCollection arguments ) {
        this.readsFiles = arguments.samFiles;
        if (arguments.downsampleFraction != null) downsamplingFraction = arguments.downsampleFraction;
        if (arguments.downsampleCoverage != null) downsampleToCoverage = arguments.downsampleCoverage;
        if (arguments.maximumReadSorts != null) maxOnFlySorts = arguments.maximumReadSorts;
        beSafe = !arguments.unsafe;
    }
}
