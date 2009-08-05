package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;

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
public class GenomeLoc implements Comparable<GenomeLoc>, Cloneable {
    private static Logger logger = Logger.getLogger(GenomeLoc.class);

    /**
     * the basic components of a genome loc, its contig index,
     * start and stop position, and (optionally) the contig name
     */
    protected final int contigIndex;
    protected final long start;
    protected final long stop;
    protected final String contigName;
    
    // --------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // --------------------------------------------------------------------------------------------------------------
    /*GenomeLoc( int contigIndex, final long start, final long stop ) {
        MAX_CONTIG = Integer.MAX_VALUE;
        if (start < 0) { throw new StingException("Bad start position " + start);}
        if (stop  < -1) { throw new StingException("Bad stop position " + stop); }    // a negative -1 indicates it's not a meaningful end position

        this.contigIndex = contigIndex;
        this.start = start;
        this.contigName = null;  // we just don't know
        this.stop = stop == -1 ? start : stop;
    }*/

    protected GenomeLoc(final SAMRecord read) {
        this(read.getHeader().getSequence(read.getReferenceIndex()).getSequenceName(), read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd());
    }

    protected GenomeLoc( final String contig, final int contigIndex, final long start, final long stop ) {
        this.contigName = contig;
        this.contigIndex = contigIndex;
        this.start = start;
        this.stop = stop;
    }

    /*GenomeLoc( final int contig, final long pos ) {
        this(contig, pos, pos );
    }
    */
    protected GenomeLoc( final GenomeLoc toCopy ) {
        this( toCopy.getContig(), toCopy.contigIndex, toCopy.getStart(), toCopy.getStop() );
    }


    /**
     * Returns true if we have a specified series of locations to process AND we are past the last
     * location in the list.  It means that, in a serial processing of the genome, that we are done.
     *
     * @param curr Current genome Location
     * @param locs a list of genomic locations
     * @return true if we are past the last location to process
     */
    public static boolean pastFinalLocation(GenomeLoc curr, List<GenomeLoc> locs) {
        return (locs.size() > 0 && curr.isPast(locs.get(locs.size() - 1)));
    }

    /**
     * A key function that returns true if the proposed GenomeLoc curr is within the list of
     * locations we are processing in this TraversalEngine
     *
     * @param curr the current location
     * @param locs a list of genomic locations
     * @return true if we should process GenomeLoc curr, otherwise false
     */
    public static boolean inLocations(GenomeLoc curr, ArrayList<GenomeLoc> locs) {
        if ( locs.size() == 0 ) {
            return true;
        } else {
            for ( GenomeLoc loc : locs ) {
                //System.out.printf("  Overlap %s vs. %s => %b%n", loc, curr, loc.overlapsP(curr));
                if (loc.overlapsP(curr))
                    return true;
            }
            return false;
        }
    }

    public static void removePastLocs(GenomeLoc curr, List<GenomeLoc> locs) {
        while ( !locs.isEmpty() && curr.isPast(locs.get(0)) ) {
            //System.out.println("At: " + curr + ", removing: " + locs.get(0));
            locs.remove(0);
        }
    }

    public static boolean overlapswithSortedLocsP(GenomeLoc curr, List<GenomeLoc> locs, boolean returnTrueIfEmpty) {
        if ( locs.isEmpty() )
            return returnTrueIfEmpty;

        // skip loci before intervals begin
        if ( curr.contigIndex < locs.get(0).contigIndex )
            return false;

        for ( GenomeLoc loc : locs ) {
            //System.out.printf("  Overlap %s vs. %s => %b%n", loc, curr, loc.overlapsP(curr));
            if ( loc.overlapsP(curr) )
                return true;
            if ( curr.compareTo(loc) < 0 )
                return false;
        }
        return false;
    }

    //
    // Accessors and setters
    //
    public final String getContig() {
        return this.contigName;
    }

    public final int getContigIndex() { return this.contigIndex; }
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

    public final boolean isUnmapped() { return this.contigIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX; }
    public final boolean throughEndOfContigP() { return this.stop == Integer.MAX_VALUE; }
    public final boolean atBeginningOfContigP() { return this.start == 1; }


    public final boolean isSingleBP() { return stop == start; }

    public final boolean disjointP(GenomeLoc that) {
        return this.contigIndex != that.contigIndex || this.start > that.stop || that.start > this.stop;
    }

    public final boolean discontinuousP(GenomeLoc that) {
        return this.contigIndex != that.contigIndex || (this.start - 1) > that.stop || (that.start - 1) > this.stop;
    }

    public final boolean overlapsP(GenomeLoc that) {
        return ! disjointP( that );
    }

    public final boolean contiguousP(GenomeLoc that) {
        return ! discontinuousP( that );
    }

    public GenomeLoc merge( GenomeLoc that ) throws StingException {
        if (!(this.contiguousP(that))) {
            throw new StingException("The two genome loc's need to be contigous");
        }

        return new GenomeLoc(getContig(), this.contigIndex,
                             Math.min(getStart(), that.getStart()),
                             Math.max( getStop(), that.getStop()) );
    }

    public final boolean containsP(GenomeLoc that) {
        return onSameContig(that) && getStart() <= that.getStart() && getStop() >= that.getStop();
    }

    public final boolean onSameContig(GenomeLoc that) {
        return (this.contigIndex == that.contigIndex);
    }

    public final int minus( final GenomeLoc that ) {
        if ( this.contigIndex == that.contigIndex )
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

    public final boolean isBefore( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == -1 || ( comparison == 0 && this.getStop() < that.getStart() ));        
    }

    public final boolean isPast( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == 1 || ( comparison == 0 && this.getStart() > that.getStop() ));
    }

    

    /**
     * Check to see whether two genomeLocs are equal.
     * Note that this implementation ignores the contigInfo object.
     * @param other Other contig to compare.
     */
    @Override
    public boolean equals(Object other) {
        if(other == null)
            return false;
        if(other instanceof GenomeLoc) {
            GenomeLoc otherGenomeLoc = (GenomeLoc)other;
            return this.contigIndex == otherGenomeLoc.contigIndex &&
                   this.start == otherGenomeLoc.start &&
                   this.stop == otherGenomeLoc.stop;
        }
        return false;
    }
    
    @Override
    public int hashCode() {
        return (int)( start << 16 + stop << 4 + contigIndex );
    }


    /**
     * Return a new GenomeLoc at this same position.
     * @return A GenomeLoc with the same contents as the current loc.
     */
    @Override
    public GenomeLoc clone() {
        return new GenomeLoc(this);
    }

    /**
     * conpare this genomeLoc's contig to another genome loc
     * @param that the genome loc to compare contigs with
     * @return 0 if equal, -1 if that.contig is greater, 1 if this.contig is greater
     */
    public final int compareContigs( GenomeLoc that ) {
        if (this.contigIndex == that.contigIndex)
            return 0;
        else if (this.contigIndex > that.contigIndex)
            return 1;
        return -1;
    }

    public int compareTo( GenomeLoc that ) {
        if ( this == that ) return 0;

        final int cmpContig = compareContigs(that);

        if ( cmpContig != 0 ) return cmpContig;
        if ( this.getStart() < that.getStart() ) return -1;
        if ( this.getStart() > that.getStart() ) return 1;

        // TODO: and error is being thrown because we are treating reads with the same start positions
        // but different stop as out of order
        //if ( this.getStop() < that.getStop() ) return -1;
        //if ( this.getStop() > that.getStop() ) return 1;
        return 0;
    }

    /**
     * How many BPs are covered by this locus?
     * @return Number of BPs covered by this locus.  According to the semantics of GenomeLoc, this should
     *         never be < 1.
     */
    public long size() {
        return stop - start + 1;
    }

}
