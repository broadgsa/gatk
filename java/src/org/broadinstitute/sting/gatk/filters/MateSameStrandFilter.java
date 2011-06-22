package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/18/11
 * Time: 4:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class MateSameStrandFilter extends ReadFilter {

    public boolean filterOut(SAMRecord read) {
        return (! read.getReadPairedFlag() ) || read.getMateUnmappedFlag() || read.getDuplicateReadFlag() ||
                read.getReadFailsVendorQualityCheckFlag() || read.getMateNegativeStrandFlag() != read.getReadNegativeStrandFlag();
    }
}
