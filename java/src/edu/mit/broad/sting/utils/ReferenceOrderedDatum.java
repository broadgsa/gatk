package edu.mit.broad.sting.utils;

import java.util.Comparator;
import java.util.HashMap;

//
// Ugly global variable defining the optional ordering of contig elements
//
/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:49:47 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReferenceOrderedDatum implements Comparable {
    public static HashMap<String, Integer> refContigOrdering = null;
    public static void setContigOrdering(HashMap<String, Integer> rco) {
        refContigOrdering = rco;
    }

    public ReferenceOrderedDatum() { }

    public abstract void parseLine(final String[] parts);

    public abstract String toString();
    public abstract String toSimpleString();
    public abstract String repl();

    public abstract String getContig();
    public abstract long getStart();
    public abstract long getStop();

    public int compareTo( Object x ) {
        if ( this == x ) return 0;

        ReferenceOrderedDatum that = (ReferenceOrderedDatum)x;

        if ( refContigOrdering != null ) {
            if ( ! refContigOrdering.containsKey(this.getContig()) ) {
                if ( ! refContigOrdering.containsKey(that.getContig()) ) {
                    // Use regular sorted order
                    int cmpContig = getContig().compareTo(that.getContig());
                    if ( cmpContig != 0 )return cmpContig;
                }
                else {
                    // this is always bigger if that is in the key set
                    return 1;   
                }
            }
            else if ( ! refContigOrdering.containsKey(that.getContig()) )
                return -1;
            else {
                assert refContigOrdering.containsKey(this.getContig()) : this;
                assert refContigOrdering.containsKey(that.getContig()) : that;

                final int thisO = refContigOrdering.get(this.getContig());
                final int thatO = refContigOrdering.get(that.getContig());
                if ( thisO < thatO ) return -1;
                if ( thisO > thatO ) return 1;
            }
        }
        else {
            int cmpContig = getContig().compareTo(that.getContig());
            if ( cmpContig != 0 )return cmpContig;
        }

        if ( this.getStart() < that.getStart() ) return -1;
        if ( this.getStart() > that.getStart() ) return 1;
        if ( this.getStop() < that.getStop() ) return -1;
        if ( this.getStop() > that.getStop() ) return 1;
        return 0;
    }
}
