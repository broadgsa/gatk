package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

/**
 * A SamRecordFilter that also depends on the header.
 */
public abstract class ReadFilter implements SamRecordFilter {
    /**
     * Sets the header for use by this filter.
     * @param engine the engine.
     */
    public void initialize(GenomeAnalysisEngine engine) {}
}
