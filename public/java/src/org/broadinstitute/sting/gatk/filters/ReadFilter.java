package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;

/**
 * A SamRecordFilter that also depends on the header.
 */
@DocumentedGATKFeature(
        groupName = "Read filters",
        summary = "GATK Engine arguments that filter or transfer incoming SAM/BAM data files" )
public abstract class ReadFilter implements SamRecordFilter {
    /**
     * Sets the header for use by this filter.
     * @param engine the engine.
     */
    public void initialize(GenomeAnalysisEngine engine) {}


    /**
     * Determines whether a pair of SAMRecord matches this filter
     *
     * @param first  the first SAMRecord to evaluate
     * @param second the second SAMRecord to evaluate
     *
     * @return true if the SAMRecords matches the filter, otherwise false
     * @throws UnsupportedOperationException when paired filter not implemented
     */
    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException("Paired filter not implemented: " + this.getClass());
    }
}
