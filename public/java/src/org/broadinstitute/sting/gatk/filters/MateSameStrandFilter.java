package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;

/**
 * Filter out reads that are not paired, have their mate unmapped, are duplicates, fail vendor quality check or both mate and read are in the same strand.
 *
 * @author chartl
 * @since 5/18/11
 */
public class MateSameStrandFilter extends ReadFilter {

    public boolean filterOut(SAMRecord read) {
        return (! read.getReadPairedFlag() ) || read.getMateUnmappedFlag() || read.getDuplicateReadFlag() ||
                read.getReadFailsVendorQualityCheckFlag() || read.getMateNegativeStrandFlag() != read.getReadNegativeStrandFlag();
    }
}
