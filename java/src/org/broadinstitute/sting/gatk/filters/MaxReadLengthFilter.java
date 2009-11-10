package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * Filters out reads whose length is >= some value.
 *
 * @author mhanna
 * @version 0.1
 */
public class MaxReadLengthFilter implements SamRecordFilter {
    @Argument(fullName = "maxReadLength", shortName = "maxRead", doc="Discard reads with length greater than the specified value", required=true)
    private int maxReadLength;
    
    public boolean filterOut(SAMRecord read) {
        // check the length
        return read.getReadLength() > maxReadLength;
    }

}
