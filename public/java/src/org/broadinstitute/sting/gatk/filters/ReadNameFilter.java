package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;

/**
 * Filter out all reads except those with this read name
 *
 * @author chartl
 * @since 9/19/11
 */
public class ReadNameFilter extends ReadFilter {
     @Argument(fullName = "readName", shortName = "rn", doc="Filter out all reads except those with this read name", required=true)
    private String readName;

    public boolean filterOut(final SAMRecord rec) {
        return ! rec.getReadName().equals(readName);
    }
}
