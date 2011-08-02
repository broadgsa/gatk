package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
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
}
