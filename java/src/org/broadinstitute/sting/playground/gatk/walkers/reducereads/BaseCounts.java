package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import org.broadinstitute.sting.utils.BaseUtils;
import com.google.java.contract.*;

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 2:55 PM
*/
final class BaseCounts {
    private final int counts[] = new int[4]; // fixme -- include - and I events

    public void incr(byte base) {
        int baseI = BaseUtils.simpleBaseToBaseIndex(base);
        if ( baseI >= 0 ) // no Ns
            counts[baseI]++;
    }

    @Ensures("BaseUtils.isRegularBase(result)")
    public byte baseWithMostCounts() {
        return BaseUtils.baseIndexToSimpleBase(maxBaseIndex());
    }

    @Ensures("result >= 0")
    public int countOfMostCommonBase() {
        return counts[maxBaseIndex()];
    }

    @Ensures("result >= 0")
    public int totalCounts() {
        int sum = 0;

        for ( int c : counts ) {
            sum += c;
        }

        return sum;
    }

    @Ensures("result >= 0 && result < counts.length")
    private int maxBaseIndex() {
        int maxI = 0;
        for ( int i = 0; i < counts.length; i++) {
            if ( counts[i] > counts[maxI] ) {
                maxI = i;
            }
        }
        return maxI;
    }

    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for ( int i = 0; i < counts.length; i++ ) {
            b.append((char)BaseUtils.baseIndexToSimpleBase(i)).append("=").append(counts[i]).append(",");
        }
        return b.toString();
    }
}
