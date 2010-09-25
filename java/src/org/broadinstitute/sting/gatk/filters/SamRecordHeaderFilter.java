package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileHeader;

/**
 * A SamRecordFilter that also depends on the header.
 */
public interface SamRecordHeaderFilter extends SamRecordFilter {
    /**
     * Sets the header for use by this filter.
     * @param header the header for use by this filter.
     */
    void setHeader(SAMFileHeader header);
}
