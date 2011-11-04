package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 9/19/11
 * Time: 4:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadNameFilter extends ReadFilter {
     @Argument(fullName = "readName", shortName = "rn", doc="Filter out all reads except those with this read name", required=true)
    private String readName;

    public boolean filterOut(final SAMRecord rec) {
        return ! rec.getReadName().equals(readName);
    }
}
