package org.broadinstitute.sting.alignment.bwa.bwt;

import org.broadinstitute.sting.utils.StingException;

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
    private Map<Byte,Integer> counts = new HashMap<Byte,Integer>();

    /**
     * Internal representation of cumulative counts, broken down by ASCII value.
     */
    private Map<Byte,Integer> cumulativeCounts = new HashMap<Byte,Integer>();

    /**
     * Create an empty Counts object with values A=0,C=0,G=0,T=0.
     */
    public Counts()
    {
        for(byte base: Bases.instance) {
            counts.put(base,0);
            cumulativeCounts.put(base,0);
        }
    }

    /**
     * Create a counts data structure with the given initial values. 
     * @param data Count data, broken down by base.
     * @param cumulative Whether the counts are cumulative, (count_G=numA+numC+numG,for example).
     */
    public Counts( int[] data, boolean cumulative ) {
        if(cumulative) {
            int priorCount = 0;
            for(byte base: Bases.instance) {
                int count = data[Bases.toPack(base)];
                counts.put(base,count-priorCount);
                cumulativeCounts.put(base,priorCount);
                priorCount = count;
            }
        }
        else {
            int priorCount = 0;
            for(byte base: Bases.instance) {
                int count = data[Bases.toPack(base)];
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
    public int[] toArray(boolean cumulative) {
        int[] countArray = new int[counts.size()];
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
            throw new StingException("Unable to clone counts object", ex);
        }
        other.counts = new HashMap<Byte,Integer>(counts);
        other.cumulativeCounts = new HashMap<Byte,Integer>(cumulativeCounts);
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
    public int get(byte base) {
        return counts.get(base);
    }

    /**
     * Gets a count of the number of bases seen before this base.
     * Note that counts in this case are cumulative.
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public int getCumulative(byte base) {
        return cumulativeCounts.get(base);
    }

    /**
     * How many total bases are represented by this count structure?
     * @return Total bases represented.
     */
    public int getTotal() {
        int accumulator = 0;
        for(byte base: Bases.instance) {
            accumulator += get(base);    
        }
        return accumulator;
    }
}
