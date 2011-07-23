package org.broadinstitute.sting.gatk.io;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;

/**
 * A writer that will allow unsorted BAM files to be written
 * and sorted on-the-fly.
 *
 * @author mhanna
 * @version 0.1
 */
public interface StingSAMFileWriter extends SAMFileWriter {
    /**
     * Writes the given custom header to SAM file output.
     * @param header The header to write.
     */
    public void writeHeader(SAMFileHeader header);

    /**
     * Set Whether the BAM file to create is actually presorted.
     * @param presorted True if the BAM file is presorted.  False otherwise.
     */    
    public void setPresorted(boolean presorted);

    /**
     * Set how many records in RAM the BAM file stores when sorting on-the-fly.
     * @param maxRecordsInRam Max number of records in RAM.
     */
    public void setMaxRecordsInRam(int maxRecordsInRam);
}