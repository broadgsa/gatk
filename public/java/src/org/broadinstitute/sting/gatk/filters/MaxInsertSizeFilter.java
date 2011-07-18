package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/2/11
 * Time: 12:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class MaxInsertSizeFilter extends ReadFilter {
    @Argument(fullName = "maxInsertSize", shortName = "maxInsert", doc="Discard reads with insert size greater than the specified value, defaults to 1000000", required=false)
    private int maxInsertSize = 1000000;

    public boolean filterOut(SAMRecord record) {
        return (record.getReadPairedFlag() && (record.getInferredInsertSize() > maxInsertSize || record.getInferredInsertSize() < -1*maxInsertSize));
    }
}
