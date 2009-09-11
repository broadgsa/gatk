package org.broadinstitute.sting.bwa;

import org.broadinstitute.sting.utils.StingException;

import java.util.EnumSet;
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
    private int[] counts = new int[EnumSet.allOf(Base.class).size()];

    /**
     * Create an empty Counts object with values A=0,C=0,G=0,T=0.
     */
    public Counts() {}

    /**
     * Create a counts data structure with the given initial values. 
     * @param data Count data, broken down by base.
     * @param cumulative Whether the counts are cumulative, (count_G=numA+numC+numG,for example).
     */
    Counts( int[] data, boolean cumulative ) {
        for( Base base: EnumSet.allOf(Base.class))
            counts[base.toPack()] = data[base.toPack()];

        // De-cumulatize data as necessary.
        if(cumulative) {
            for( int i = EnumSet.allOf(Base.class).size()-1; i > 0; i-- )
                counts[i] -= counts[i-1];
        }
    }

    /**
     * Convert to an array for persistence.
     * @param cumulative Use a cumulative representation.
     * @return Array of count values.
     */
    public int[] toArray(boolean cumulative) {
        int[] countArray = counts.clone();
        if(cumulative) {
            for( int i = 1; i < counts.length; i++ )
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
        other.counts = new int[counts.length];
        System.arraycopy(counts,0,other.counts,0,counts.length);
        return other;
    }

    /**
     * Increment the number of bases seen at the given location.
     * @param base Base to increment.
     */
    public void increment(Base base) {
        counts[base.toPack()]++;        
    }

    /**
     * Gets a count of the number of bases seen at a given location.
     * Note that counts in this case are not cumulative (counts for A,C,G,T
     * are independent).
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public int get(Base base) {
        return counts[base.toPack()];
    }

    /**
     * Gets a count of the number of bases of each seen type.
     * Note that counts in this case are cumulative (counts for A,C,G,T
     * are independent).
     * @param base Base for which to query counts.
     * @return Number of bases of this type seen.
     */
    public int getCumulative(Base base) {
        int accum = 0;
        for(int i = 0; i <= base.toPack(); i++)
            accum += counts[i];
        return accum;
    }

    /**
     * How many total bases are represented by this count structure?
     * @return Total bases represented.
     */
    public int getTotal() {
        int accumulator = 0;
        for( int count : counts )
            accumulator += count;
        return accumulator;
    }
}
