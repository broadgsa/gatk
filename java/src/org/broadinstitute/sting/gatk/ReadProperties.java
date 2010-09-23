package org.broadinstitute.sting.gatk;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
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
public class ReadProperties {
    private List<SAMReaderID> readers = null;
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.STRICT;
    private Integer readBufferSize = null;
    private DownsamplingMethod downsamplingMethod = null;
    private ValidationExclusion exclusionList = null;
    private Collection<SamRecordFilter> supplementalFilters = null;
    private boolean includeReadsWithDeletionAtLoci = false;
    private boolean generateExtendedEvents = false; // do we want to generate additional piles of "extended" events (indels)
// immediately after the reference base such event is associated with?


    /**
     * Return true if the walker wants to see reads that contain deletions when looking at locus pileups
     * 
     * @return
     */
    public boolean includeReadsWithDeletionAtLoci() {
        return includeReadsWithDeletionAtLoci;
    }

    /**
     * Return true if the walker wants to see additional piles of "extended" events (indels). An indel is associated,
     * by convention, with the reference base immediately preceding the insertion/deletion, and if this flag is set
     * to 'true', any locus with an indel associated with it will cause exactly two subsequent calls to walker's map(): first call
     * will be made with a "conventional" base pileup, the next call will be made with a pileup of extended (indel/noevent)
     * events.
     * @return
     */
    public boolean generateExtendedEvents() {
        return generateExtendedEvents;
    }

    /**
     * Gets a list of the files acting as sources of reads.
     * @return A list of files storing reads data.
     */
    public List<SAMReaderID> getSAMReaderIDs() {
        return readers;
    }

    /**
     * How strict should validation be?
     * @return Stringency of validation.
     */
    public SAMFileReader.ValidationStringency getValidationStringency() {
        return validationStringency;
    }

    /**
     * Gets a list of the total number of reads that the sharding system should buffer per BAM file.
     * @return
     */
    public Integer getReadBufferSize() {
        return readBufferSize;
    }

    /**
     * Gets the method and parameters used when downsampling reads.
     * @return Downsample fraction.
     */
    public DownsamplingMethod getDownsamplingMethod() {
        return downsamplingMethod;
    }

    /**
     * Return whether to 'verify' the reads as we pass through them.
     * @return Whether to verify the reads.
     */
    public ValidationExclusion getValidationExclusionList() {
        return exclusionList;
    }

    public Collection<SamRecordFilter> getSupplementalFilters() {
        return supplementalFilters;
    }

    /**
     * Simple constructor for unit testing.
     * @param readsFiles List of reads files to open.
     */
    public ReadProperties( List<SAMReaderID> readsFiles ) {
        this.readers = readsFiles;
        this.downsamplingMethod = DownsamplingMethod.NONE;
        this.supplementalFilters = new ArrayList<SamRecordFilter>();
        this.exclusionList = new ValidationExclusion();
    }

    /**
     * Extract the command-line arguments having to do with reads input
     * files and store them in an easy-to-work-with package.  Constructor
     * is package protected.
     * @param samFiles list of reads files.
     * @param strictness Stringency of reads file parsing.
     * @param readBufferSize Number of reads to hold in memory per BAM.
     * @param downsamplingMethod Method for downsampling reads at a given locus.
     * @param exclusionList what safety checks we're willing to let slide
     * @param supplementalFilters additional filters to dynamically apply.
     * @param generateExtendedEvents if true, the engine will issue an extra call to walker's map() with
     *        a pile of indel/noevent extended events at every locus with at least one indel associated with it
     *        (in addition to a "regular" call to map() at this locus performed with base pileup)
     * @param includeReadsWithDeletionAtLoci if 'true', the base pileups sent to the walker's map() method
     *         will explicitly list reads with deletion over the current reference base; otherwise, only observed
     *        bases will be seen in the pileups, and the deletions will be skipped silently.
     */
    public ReadProperties( List<SAMReaderID> samFiles,
           SAMFileReader.ValidationStringency strictness,
           Integer readBufferSize,
           DownsamplingMethod downsamplingMethod,
           ValidationExclusion exclusionList,
           Collection<SamRecordFilter> supplementalFilters,
           boolean includeReadsWithDeletionAtLoci,
           boolean generateExtendedEvents) {
        this.readers = samFiles;
        this.readBufferSize = readBufferSize;
        this.validationStringency = strictness;
        this.downsamplingMethod = downsamplingMethod;
        this.exclusionList = exclusionList == null ? new ValidationExclusion() : exclusionList;
        this.supplementalFilters = supplementalFilters;
        this.includeReadsWithDeletionAtLoci = includeReadsWithDeletionAtLoci;
        this.generateExtendedEvents = generateExtendedEvents;
    }
}
