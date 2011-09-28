package org.broadinstitute.sting.utils.sam;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMRecord;

import java.util.Comparator;

public class AlignmentStartWithNoTiesComparator implements Comparator<SAMRecord> {
    @Requires("c1 >= 0 && c2 >= 0")
    @Ensures("result == 0 || result == 1 || result == -1")
    private int compareContigs(int c1, int c2) {
        if (c1 == c2)
            return 0;
        else if (c1 > c2)
            return 1;
        return -1;
    }

    @Requires("r1 != null && r2 != null")
    @Ensures("result == 0 || result == 1 || result == -1")
    public int compare(SAMRecord r1, SAMRecord r2) {
        int result;

        if (r1 == r2)
            result = 0;

        else if (r1.getReadUnmappedFlag())
            result = 1;
        else if (r2.getReadUnmappedFlag())
            result = -1;
        else {
            final int cmpContig = compareContigs(r1.getReferenceIndex(), r2.getReferenceIndex());

            if (cmpContig != 0)
                result = cmpContig;

            else {
                if (r1.getAlignmentStart() < r2.getAlignmentStart()) result = -1;
                else result = 1;
            }
        }

        return result;
    }
}