package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import org.broadinstitute.sting.utils.BaseUtils;
import com.google.java.contract.*;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 4/8/11
* Time: 2:55 PM
*/
final class BaseCounts {
    private final int counts[] = new int[4]; // todo -- fixme -- include - and I events

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
