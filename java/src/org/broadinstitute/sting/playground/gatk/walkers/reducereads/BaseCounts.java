package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import org.broadinstitute.sting.alignment.reference.bwt.Counts;
import org.broadinstitute.sting.utils.BaseUtils;
import com.google.java.contract.*;

import java.util.EnumMap;
import java.util.Map;

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
    public final static BaseIndex MAX_BASE_INDEX_WITH_NO_COUNTS = BaseIndex.A;
    public final static byte MAX_BASE_WITH_NO_COUNTS = MAX_BASE_INDEX_WITH_NO_COUNTS.getByte();

    private final Map<BaseIndex, Integer> counts = new EnumMap<BaseIndex, Integer>(BaseIndex.class); // todo -- fixme -- include - and I events
    {
        for ( BaseIndex i : BaseIndex.values() )
            counts.put(i,0);
    }

    @Ensures("totalCount() == old(totalCount()) || totalCount() == old(totalCount()) + 1")
    public void incr(byte base) {
        BaseIndex i = BaseIndex.byteToBase(base);
        if ( i != null ) // no Ns
            counts.put(i, counts.get(i) + 1);
    }

    public byte baseWithMostCounts() {
        return maxBaseIndex().getByte();
    }

    @Ensures("result >= 0")
    public int countOfMostCommonBase() {
        return counts.get(maxBaseIndex());
    }

    @Ensures("result >= 0")
    public int totalCount() {
        int sum = 0;

        for ( int c : counts.values() ) {
            sum += c;
        }

        return sum;
    }

    @Ensures({
            "result != null",
            "totalCount() != 0 || result == MAX_BASE_INDEX_WITH_NO_COUNTS"})
    private BaseIndex maxBaseIndex() {
        BaseIndex maxI = MAX_BASE_INDEX_WITH_NO_COUNTS;
        for ( BaseIndex i : counts.keySet() ) {
            if ( counts.get(i) > counts.get(maxI) ) {
                maxI = i;
            }
        }
        return maxI;
    }

    @Ensures("result != null")
    public String toString() {
        StringBuilder b = new StringBuilder();
        for ( Map.Entry<BaseIndex,Integer> elt : counts.entrySet() ) {
            b.append(elt.toString()).append("=").append(elt.getValue()).append(",");
        }
        return b.toString();
    }

    private enum BaseIndex {
        A ( 'A', 0 ),
        C ( 'C', 1 ),
        G ( 'G', 2 ),
        T ( 'T', 3 ),
        D ( 'D', 4 ),
        I ( 'I', 5 ); // insertion to the right of the base

        final byte b;
        final int index;
        private BaseIndex(char base, int index) {
            this.b = (byte)base;
            this.index = index;
        }

        public byte getByte() { return b; }

        public static final BaseIndex byteToBase(final byte base) {
            switch (base) {
                case 'A': return A;
                case 'C': return C;
                case 'G': return G;
                case 'T': return T;
                case 'D': return D;
                case 'I': return I;
                default: return null;
            }
        }
    }
}
