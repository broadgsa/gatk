/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.HashMap;
import java.util.Map;

/**
 * Counts of how many bases of each type have been seen.
 *
 * @author mhanna
 * @version 0.1
 */
public class Counts implements Cloneable {
    /**
     * Internal representation of counts, broken down by ASCII value.
     */
    private Map<Byte,Long> counts = new HashMap<Byte,Long>();

    /**
     * Internal representation of cumulative counts, broken down by ASCII value.
     */
    private Map<Byte,Long> cumulativeCounts = new HashMap<Byte,Long>();

    /**
     * Create an empty Counts object with values A=0,C=0,G=0,T=0.
     */
    public Counts()
    {
        for(byte base: Bases.instance) {
            counts.put(base,0L);
            cumulativeCounts.put(base,0L);
        }
    }

    /**
     * Create a counts data structure with the given initial values. 
     * @param data Count data, broken down by base.
     * @param cumulative Whether the counts are cumulative, (count_G=numA+numC+numG,for example).
     */
    public Counts( long[] data, boolean cumulative ) {
        if(cumulative) {
            long priorCount = 0;
            for(byte base: Bases.instance) {
                long count = data[Bases.toPack(base)];
                counts.put(base,count-priorCount);
                cumulativeCounts.put(base,priorCount);
                priorCount = count;
            }
        }
        else {
            long priorCount = 0;
            for(byte base: Bases.instance) {
                long count = data[Bases.toPack(base)];
                counts.put(base,count);
                cumulativeCounts.put(base,priorCount);
                priorCount += count;
            }
        }
    }

    /**
     * Convert to an array for persistence.
     * @param cumulative Use a cumulative representation.
     * @return Array of count values.
     */
    public long[] toArray(boolean cumulative) {
        long[] countArray = new long[counts.size()];
        if(cumulative) {
            int index = 0;
            boolean first = true;
            for(byte base: Bases.instance) {
                if(first) {
                    first = false;
                    continue;
                }
                countArray[index++] = getCumulative(base);
            }
            countArray[countArray.length-1] = getTotal();
        }
        else {
            int index = 0;
            for(byte base: Bases.instance)
                countArray[index++] = counts.get(base);
        }
        return countArray;
    }

    /**
     * Create a unique copy of the current object.
     * @return A duplicate of this object.
     */
    public Counts clone() {
        Counts other;
        try {
            other = (Counts)super.clone();
        }
        catch(CloneNotSupportedException ex) {
            throw new ReviewedGATKException("Unable to clone counts object", ex);
        }
        other.counts = new HashMap<Byte,Long>(counts);
        other.cumulativeCounts = new HashMap<Byte,Long>(cumulativeCounts);
        return other;
    }

    /**
     * Increment the number of bases seen at the given location.
     * @param base Base to increment.
     */
    public void increment(byte base) {
        counts.put(base,counts.get(base)+1);
        boolean increment = false;
        for(byte cumulative: Bases.instance) {
            if(increment) cumulativeCounts.put(cumulative,cumulativeCounts.get(cumulative)+1);
            increment |= (cumulative == base);
        }
    }

    /**
     * Gets a count of the number of bases seen at a given location.
     * Note that counts in this case are not cumulative (counts for A,C,G,T
     * are independent).
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public long get(byte base) {
        return counts.get(base);
    }

    /**
     * Gets a count of the number of bases seen before this base.
     * Note that counts in this case are cumulative.
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public long getCumulative(byte base) {
        return cumulativeCounts.get(base);
    }

    /**
     * How many total bases are represented by this count structure?
     * @return Total bases represented.
     */
    public long getTotal() {
        int accumulator = 0;
        for(byte base: Bases.instance) {
            accumulator += get(base);    
        }
        return accumulator;
    }
}
