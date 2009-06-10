package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import net.sf.samtools.SAMFileReader;
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
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.STRICT;
    private Double downsamplingFraction = null;
    private Integer downsampleToCoverage = null;
    private Integer maxOnFlySorts = null;
    private Boolean beSafe = null;
    private Boolean filterZeroMappingQualityReads = null;

    /**
     * Gets a list of the files acting as sources of reads.
     * @return A list of files storing reads data.
     */
    public List<File> getReadsFiles() {
        return readsFiles;
    }

    /**
     * How strict should validation be?
     * @return Stringency of validation.
     */
    public SAMFileReader.ValidationStringency getValidationStringency() {
        return validationStringency;
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

    public Boolean getFilterZeroMappingQualityReads() {
        return filterZeroMappingQualityReads;
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
     * @param samFiles list of reads files.
     * @param strictness Stringency of reads file parsing.
     * @param downsampleFraction fraction of reads to downsample.
     * @param downsampleCoverage downsampling per-locus.
     * @param maxOnFlySorts how many sorts to perform on-the-fly.
     * @param beSafe Whether to enable safety checking.
     * @param filterZeroMappingQualityReads whether to filter zero mapping quality reads.
     */
    Reads( List<File> samFiles,
           SAMFileReader.ValidationStringency strictness,
           Double downsampleFraction,
           Integer downsampleCoverage,
           Integer maxOnFlySorts,
           Boolean beSafe,
           Boolean filterZeroMappingQualityReads ) {
        this.readsFiles = samFiles;
        this.validationStringency = strictness;
        this.downsamplingFraction = downsampleFraction;
        this.downsampleToCoverage = downsampleCoverage;
        this.maxOnFlySorts = maxOnFlySorts;
        this.beSafe = beSafe;
        this.filterZeroMappingQualityReads = filterZeroMappingQualityReads;
    }
}
