package org.broadinstitute.sting.gatk;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.baq.BAQ;

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
    private Collection<SAMReaderID> readers = null;
    private SAMFileHeader header = null;
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.STRICT;
    private Integer readBufferSize = null;
    private DownsamplingMethod downsamplingMethod = null;
    private ValidationExclusion exclusionList = null;
    private Collection<ReadFilter> supplementalFilters = null;
    private boolean includeReadsWithDeletionAtLoci = false;
    private boolean useOriginalBaseQualities = false;
    private boolean generateExtendedEvents = false;
    private BAQ.CalculationMode cmode = BAQ.CalculationMode.OFF;
    private BAQ.QualityMode qmode = BAQ.QualityMode.DONT_MODIFY;
    IndexedFastaSequenceFile refReader = null; // read for BAQ, if desired
    private byte defaultBaseQualities;

    // do we want to generate additional piles of "extended" events (indels)
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
    public Collection<SAMReaderID> getSAMReaderIDs() {
        return readers;
    }

    /**
     * Gets the sam file header
     * @return the sam file header
     */
    public SAMFileHeader getHeader() {
        return header;
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

    public Collection<ReadFilter> getSupplementalFilters() {
        return supplementalFilters;
    }

    /**
     * Return whether to use original base qualities.
     * @return Whether to use original base qualities.
     */
    public boolean useOriginalBaseQualities() {
        return useOriginalBaseQualities;
    }


    public BAQ.QualityMode getBAQQualityMode() { return qmode; }
    public BAQ.CalculationMode getBAQCalculationMode() { return cmode; }

    public IndexedFastaSequenceFile getRefReader() {
        return refReader;
    }

    /**
     * @return Default base quality value to fill reads missing base quality information.
     */
    public byte defaultBaseQualities() {
        return defaultBaseQualities;
    }

    /**
     * Extract the command-line arguments having to do with reads input
     * files and store them in an easy-to-work-with package.  Constructor
     * is package protected.
     * @param samFiles list of reads files.
     * @param header sam file header.
     * @param useOriginalBaseQualities True if original base qualities should be used.
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
     * @param cmode How should we apply the BAQ calculation to the reads?
     * @param qmode How should we apply the BAQ calculation to the reads?
     * @param refReader if applyBAQ is true, must be a valid pointer to a indexed fasta file reads so we can get the ref bases for BAQ calculation
     * @param defaultBaseQualities if the reads have incomplete quality scores, set them all to defaultBaseQuality.
     */
    public ReadProperties( Collection<SAMReaderID> samFiles,
           SAMFileHeader header,
           boolean useOriginalBaseQualities,
           SAMFileReader.ValidationStringency strictness,
           Integer readBufferSize,
           DownsamplingMethod downsamplingMethod,
           ValidationExclusion exclusionList,
           Collection<ReadFilter> supplementalFilters,
           boolean includeReadsWithDeletionAtLoci,
           boolean generateExtendedEvents,
           BAQ.CalculationMode cmode,
           BAQ.QualityMode qmode,
           IndexedFastaSequenceFile refReader,
           byte defaultBaseQualities) {
        this.readers = samFiles;
        this.header = header;
        this.readBufferSize = readBufferSize;
        this.validationStringency = strictness;
        this.downsamplingMethod = downsamplingMethod == null ? DownsamplingMethod.NONE : downsamplingMethod;
        this.exclusionList = exclusionList == null ? new ValidationExclusion() : exclusionList;
        this.supplementalFilters = supplementalFilters;
        this.includeReadsWithDeletionAtLoci = includeReadsWithDeletionAtLoci;
        this.generateExtendedEvents = generateExtendedEvents;
        this.useOriginalBaseQualities = useOriginalBaseQualities;
        this.cmode = cmode;
        this.qmode = qmode;
        this.refReader = refReader;
        this.defaultBaseQualities = defaultBaseQualities;
    }
}
