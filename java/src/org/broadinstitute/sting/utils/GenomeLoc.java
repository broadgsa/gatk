package org.broadinstitute.sting.utils;

import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 2, 2009
 * Time: 8:50:11 AM
 *
 * Genome location representation.  It is *** 1 *** based
 *
 *
 */
public class GenomeLoc implements Comparable<GenomeLoc> {
    private String contig;
    private long start;
    private long stop;

    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    public static Map<String, Integer> refContigOrdering = null;
    public static HashMap<String, String> interns = null;

    public static void setContigOrdering(Map<String, Integer> rco) {
        refContigOrdering = rco;
        interns = new HashMap<String, String>();
        for ( String contig : rco.keySet() )
            interns.put( contig, contig );
    }

    public GenomeLoc( String contig, final long start, final long stop ) {
        if ( interns != null )
            contig = interns.get(contig);

        this.contig = contig;
        this.start = start;
        this.stop = stop;
    }

    public GenomeLoc( final String contig, final long pos ) {
        this( contig, pos, pos );
    }

    public GenomeLoc( final GenomeLoc toCopy ) {
        this( new String(toCopy.getContig()), toCopy.getStart(), toCopy.getStop() );
    }

    //
    // Parsing string representations
    //
    private static long parsePosition( final String pos ) {
        String x = pos.replaceAll(",", "");
        return Long.parseLong(x);
    }

    public static GenomeLoc parseGenomeLoc( final String str ) {
        // Ôchr2Õ, Ôchr2:1000000Õ or Ôchr2:1,000,000-2,000,000Õ
        System.out.printf("Parsing location '%s'%n", str);

        final Pattern regex1 = Pattern.compile("([\\w&&[^:]]+)$");             // matches case 1
        final Pattern regex2 = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)$");      // matches case 2
        final Pattern regex3 = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)-([\\d,]+)$");// matches case 3

        String contig = null;
        long start = 1;
        long stop = Integer.MAX_VALUE;
        boolean bad = false;

        Matcher match1 = regex1.matcher(str);
        Matcher match2 = regex2.matcher(str);
        Matcher match3 = regex3.matcher(str);

        try {
            if ( match1.matches() ) {
                contig = match1.group(1);
            }
            else if ( match2.matches() ) {
                contig = match2.group(1);
                start = parsePosition(match2.group(2));
            }
            else if ( match3.matches() ) {
                contig = match3.group(1);
                start = parsePosition(match3.group(2));
                stop = parsePosition(match3.group(3));

                if ( start > stop )
                    bad = true;
            }
            else {
                bad = true;
            }
        } catch ( Exception e ) {
            bad = true;
        }

        if ( bad ) {
            throw new RuntimeException("Invalid Genome Location string: " + str);
        }

        GenomeLoc loc = new GenomeLoc(contig, start, stop);
        System.out.printf("  => Parsed location '%s' into %s%n", str, loc);

        return loc;
    }

    //
    // Accessors and setters
    //
    public final String getContig() { return this.contig; }
    public final long getStart()    { return this.start; }
    public final long getStop()     { return this.stop; }
    public final String toString()  {
        if ( throughEndOfContigP() && atBeginningOfContigP() )
            return getContig();
        else if ( throughEndOfContigP() || getStart() == getStop() )
            return String.format("%s:%d", getContig(), getStart());
        else
            return String.format("%s:%d-%d", getContig(), getStart(), getStop());
    }


    public final boolean isUnmapped() { return this.contig == null; }
    public final boolean throughEndOfContigP() { return this.stop == Integer.MAX_VALUE; }
    public final boolean atBeginningOfContigP() { return this.start == 1; }

    public void setContig(String contig) {
        this.contig = contig;
    }

    public void setStart(long start) {
        this.start = start;
    }
    public void setStop(long stop) {
        this.stop = stop;
    }

    public final boolean isSingleBP() { return stop == start; }

    public final boolean disjointP(GenomeLoc that) {
        if ( compareContigs(this.contig, that.contig) != 0 ) return true;   // different chromosomes
        if ( this.start > that.stop ) return true;                          // this guy is past that
        if ( that.start > this.stop ) return true;                          // that guy is past our start
        return false;
    }

    public final boolean overlapsP(GenomeLoc that) {
        return ! disjointP( that );
    }

    public final boolean onSameContig(GenomeLoc that) {
        return this.contig.equals(that.contig);
    }

    public final int minus( final GenomeLoc that ) {
        if ( this.getContig().equals(that.getContig()) )
            return (int) (this.getStart() - that.getStart());
        else
            return Integer.MAX_VALUE;
    }

    public final int distance( final GenomeLoc that ) {
        return Math.abs(minus(that));
    }    

    public final boolean isBetween( final GenomeLoc left, final GenomeLoc right ) {
        return this.compareTo(left) > -1 && this.compareTo(right) < 1;
    }

    public final void incPos() {
        incPos(1);
    }
    public final void incPos(long by) {
        this.start += by;
        this.stop += by;
    }

    public final GenomeLoc nextLoc() {
        GenomeLoc n = new GenomeLoc(this);
        n.incPos();
        return n;
    }
    //
    // Comparison operations
    //
    public static int compareContigs( final String thisContig, final String thatContig ) 
    {
        if ( thisContig == thatContig )
        {
            // Optimization.  If the pointers are equal, then the contigs are equal.
            return 0;
        }

        assert refContigOrdering.containsKey(thisContig);// : this;
        assert refContigOrdering.containsKey(thatContig);// : that;

        if ( refContigOrdering != null ) 
        {
            if ( ! refContigOrdering.containsKey(thisContig) ) 
            {
                if ( ! refContigOrdering.containsKey(thatContig) ) 
                {
                    // Use regular sorted order
                    return thisContig.compareTo(thatContig);
                }
                else 
                {
                    // this is always bigger if that is in the key set
                    return 1;
                }
            }
            else if ( ! refContigOrdering.containsKey(thatContig) )
            {
                return -1;
            }
            else 
            {
                assert refContigOrdering.containsKey(thisContig);// : this;
                assert refContigOrdering.containsKey(thatContig);// : that;

                final int thisO = refContigOrdering.get(thisContig);
                final int thatO = refContigOrdering.get(thatContig);

                if ( thisO < thatO ) return -1;
                if ( thisO > thatO ) return 1;
                return 0;
            }
        }
        else 
        {
            return thisContig.compareTo(thatContig);
        }
    }

    public int compareContigs( GenomeLoc that ) {
        return compareContigs( this.contig, that.contig );
    }


    public int compareTo( GenomeLoc that ) {
        if ( this == that ) return 0;

        final int cmpContig = compareContigs( this.getContig(), that.getContig() );
        if ( cmpContig != 0 ) return cmpContig;
        if ( this.getStart() < that.getStart() ) return -1;
        if ( this.getStart() > that.getStart() ) return 1;
        if ( this.getStop() < that.getStop() ) return -1;
        if ( this.getStop() > that.getStop() ) return 1;
        return 0;
    }
}
