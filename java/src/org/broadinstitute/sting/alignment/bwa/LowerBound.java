package org.broadinstitute.sting.alignment.bwa;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.alignment.bwa.bwt.Base;
import org.broadinstitute.sting.alignment.bwa.bwt.BWT;

/**
 * At any point along the given read, what is a good lower bound for the
 * total number of differences?
 *
 * @author mhanna
 * @version 0.1
 */
public class LowerBound {
    /**
     * Lower bound of the suffix array.
     */
    public final int loIndex;

    /**
     * Upper bound of the suffix array.
     */
    public final int hiIndex;

    /**
     * Width of the bwt from loIndex -> hiIndex, inclusive.
     */
    public final int width;

    /**
     * The lower bound at the given point.
     */
    public final int value;

    /**
     * Create a new lower bound with the given value.
     * @param value Value for the lower bound at this site.
     */
    private LowerBound(int loIndex, int hiIndex, int value) {
        this.loIndex = loIndex;
        this.hiIndex = hiIndex;
        this.width = hiIndex - loIndex + 1;
        this.value = value;
    }

    /**
     * Create a non-optimal bound according to the algorithm specified in Figure 3 of the BWA paper.
     */
    public static List<LowerBound> create( byte[] bases, BWT bwt ) {
        List<LowerBound> bounds = new ArrayList<LowerBound>();

        int loIndex = 0, hiIndex = bwt.length(), mismatches = 0;
        for( int i = bases.length-1; i >= 0; i-- ) {
            Base base = Base.fromASCII(bases[i]);
            loIndex = bwt.counts(base) + bwt.occurrences(base,loIndex-1) + 1;
            hiIndex = bwt.counts(base) + bwt.occurrences(base,hiIndex);
            if( loIndex > hiIndex ) {
                loIndex = 0;
                hiIndex = bwt.length();
                mismatches++;
            }
            bounds.add(0,new LowerBound(loIndex,hiIndex,mismatches));
        }

        return bounds;
    }

    /**
     * Create a string representation of this bound.
     * @return String version of this bound.
     */
    public String toString() {
        return String.format("LowerBound: w = %d, value = %d",width,value);
    }
}
