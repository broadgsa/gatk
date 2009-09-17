package org.broadinstitute.sting.bwa;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

/**
 * At any point along the given read, what is a good lower bound for the
 * total number of differences?
 *
 * @author mhanna
 * @version 0.1
 */
public class LowerBound {
    /**
     * The lower bound at the given point.
     */
    public final int value;

    /**
     * Create a new lower bound with the given value.
     * @param value Value for the lower bound at this site.
     */
    private LowerBound(int value) {
        this.value = value;
    }

    /**
     * Create a non-optimal bound according to the algorithm specified in Figure 3 of the BWA paper.
     */
    public static List<LowerBound> create( SAMRecord read, BWT reverseBWT ) {
        List<LowerBound> bounds = new ArrayList<LowerBound>();

        int loIndex = 0, hiIndex = reverseBWT.length(), mismatches = 0;
        for( int i = 0; i < read.getReadBases().length; i++ ) {
            Base base = Base.fromASCII(read.getReadBases()[i]);
            loIndex = reverseBWT.counts(base) + reverseBWT.occurrences(base,loIndex-1) + 1;
            hiIndex = reverseBWT.counts(base) + reverseBWT.occurrences(base,hiIndex);
            if( loIndex > hiIndex ) {
                loIndex = 0;
                hiIndex = reverseBWT.length();
                mismatches++;
            }
            bounds.add(new LowerBound(mismatches));
        }

        return bounds;
    }
}
