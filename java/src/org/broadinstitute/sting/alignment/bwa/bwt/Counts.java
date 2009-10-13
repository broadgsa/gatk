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
     * Internal representation of counts, broken down by pack value.
     */
    private Map<Byte,Integer> counts = new HashMap<Byte,Integer>();

    /**
     * Create an empty Counts object with values A=0,C=0,G=0,T=0.
     */
    public Counts() {}

    /**
     * Create a counts data structure with the given initial values. 
     * @param data Count data, broken down by base.
     * @param cumulative Whether the counts are cumulative, (count_G=numA+numC+numG,for example).
     */
    public Counts( int[] data, boolean cumulative ) {
        for( byte base: Bases.instance)
            counts.put(base,data[Bases.toPack(base)]);

        // De-cumulatize data as necessary.
        if(cumulative) {
            int previousCount = 0;
            for( byte base: Bases.instance ) {
                int count = counts.get(base);
                counts.put(base,count-previousCount);
                previousCount = count;
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
        for(byte base: Bases.instance)
            countArray[Bases.toPack(base)] = counts.get(base);
        if(cumulative) {
            for( int i = 1; i < countArray.length; i++ )
                countArray[i] += countArray[i-1];
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
        return other;
    }

    /**
     * Increment the number of bases seen at the given location.
     * @param base Base to increment.
     */
    public void increment(byte base) {
        counts.put(base,counts.get(base)+1);
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
     * Gets a count of the number of bases of each seen type.
     * Note that counts in this case are cumulative (counts for A,C,G,T
     * are independent).
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public int getCumulative(byte base) {
        int accum = 0;
        for( byte current: Bases.allOf() ) {
            if(base == current) break;
            accum += counts.get(current);
        }
        return accum;
    }

    /**
     * How many total bases are represented by this count structure?
     * @return Total bases represented.
     */
    public int getTotal() {
        int accumulator = 0;
        for( int count : counts.values() )
            accumulator += count;
        return accumulator;
    }
}
