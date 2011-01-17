package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.List;
import java.io.Serializable;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 2, 2009
 * Time: 8:50:11 AM
 *
 * Genome location representation.  It is *** 1 *** based closed.
 *
 *
 */
public class GenomeLoc implements Comparable<GenomeLoc>, Cloneable, Serializable {
    /**
     * the basic components of a genome loc, its contig index,
     * start and stop position, and (optionally) the contig name
     */
    protected final int contigIndex;
    protected final int start;
    protected final int stop;
    protected final String contigName;

    /**
     * A static constant to use when referring to the unmapped section of a datafile
     * file.  The unmapped region cannot be subdivided.  Only this instance of
     * the object may be used to refer to the region, as '==' comparisons are used
     * in comparators, etc.
     */
    // TODO - WARNING WARNING WARNING code somehow depends on the name of the contig being null!
    public static final GenomeLoc UNMAPPED = new GenomeLoc(null,-1,0,0);
    public static final boolean isUnmapped(GenomeLoc loc) {
        return loc == UNMAPPED;
    }
    
    // --------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // --------------------------------------------------------------------------------------------------------------

    protected GenomeLoc(final SAMRecord read) {
        this(read.getHeader().getSequence(read.getReferenceIndex()).getSequenceName(), read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd());
    }

    protected GenomeLoc( final String contig, final int contigIndex, final int start, final int stop ) {
        this.contigName = contig;
        this.contigIndex = contigIndex;
        this.start = start;
        this.stop = stop;
    }

    /**
     * Return a new GenomeLoc at this same position.
     * @return A GenomeLoc with the same contents as the current loc.
     */
    @Override
    public GenomeLoc clone() {
        return new GenomeLoc(getContig(),getContigIndex(),getStart(),getStop());
    }

    //
    // Accessors and setters
    //
    public final String getContig() {
        return this.contigName;
    }

    public final int getContigIndex() { return this.contigIndex; }
    public final int getStart()    { return this.start; }
    public final int getStop()     { return this.stop; }
    public final String toString()  {
        if(GenomeLoc.isUnmapped(this)) return "unmapped";
        if ( throughEndOfContigP() && atBeginningOfContigP() )
            return getContig();
        else if ( throughEndOfContigP() || getStart() == getStop() )
            return String.format("%s:%d", getContig(), getStart());
        else
            return String.format("%s:%d-%d", getContig(), getStart(), getStop());
    }
    private boolean throughEndOfContigP() { return this.stop == Integer.MAX_VALUE; }
    private boolean atBeginningOfContigP() { return this.start == 1; }    

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

    public GenomeLoc merge( GenomeLoc that ) throws ReviewedStingException {
        if(GenomeLoc.isUnmapped(this) || GenomeLoc.isUnmapped(that)) {
            if(! GenomeLoc.isUnmapped(this) || !GenomeLoc.isUnmapped(that))
                throw new ReviewedStingException("Tried to merge a mapped and an unmapped genome loc");
            return UNMAPPED;
        }

        if (!(this.contiguousP(that))) {
            throw new ReviewedStingException("The two genome loc's need to be contigous");
        }

        return new GenomeLoc(getContig(), this.contigIndex,
                             Math.min(getStart(), that.getStart()),
                             Math.max( getStop(), that.getStop()) );
    }

    public GenomeLoc intersect( GenomeLoc that ) throws ReviewedStingException {
        if(GenomeLoc.isUnmapped(this) || GenomeLoc.isUnmapped(that)) {
            if(! GenomeLoc.isUnmapped(this) || !GenomeLoc.isUnmapped(that))
                throw new ReviewedStingException("Tried to intersect a mapped and an unmapped genome loc");
            return UNMAPPED;
        }

        if (!(this.overlapsP(that))) {
            throw new ReviewedStingException("GenomeLoc::intersect(): The two genome loc's need to overlap");
        }

        return new GenomeLoc(getContig(), this.contigIndex,
                             Math.max(getStart(), that.getStart()),
                             Math.min( getStop(), that.getStop()) );
    }

    public final boolean containsP(GenomeLoc that) {
        return onSameContig(that) && getStart() <= that.getStart() && getStop() >= that.getStop();
    }

    public final boolean onSameContig(GenomeLoc that) {
        return (this.contigIndex == that.contigIndex);
    }

    public final int minus( final GenomeLoc that ) {
        if ( this.onSameContig(that) )
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

    /**
     * Tests whether this contig is completely before contig 'that'.
     * @param that Contig to test against.
     * @return true if this contig ends before 'that' starts; false if this is completely after or overlaps 'that'.
     */
    public final boolean isBefore( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == -1 || ( comparison == 0 && this.getStop() < that.getStart() ));        
    }

    /**
     * Tests whether this contig is completely after contig 'that'.
     * @param that Contig to test against.
     * @return true if this contig starts after 'that' ends; false if this is completely before or overlaps 'that'.
     */
    public final boolean isPast( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == 1 || ( comparison == 0 && this.getStart() > that.getStop() ));
    }

    // Return the minimum distance between any pair of bases in this and that GenomeLocs:
    public final int minDistance( final GenomeLoc that ) {
        if (!this.onSameContig(that))
            return Integer.MAX_VALUE;

        int minDistance;
        if (this.isBefore(that))
            minDistance = distanceFirstStopToSecondStart(this, that);
        else if (that.isBefore(this))
            minDistance = distanceFirstStopToSecondStart(that, this);
        else // this and that overlap [and possibly one contains the other]:
            minDistance = 0;

        return minDistance;
    }

    private static int distanceFirstStopToSecondStart(GenomeLoc locFirst, GenomeLoc locSecond) {
        return (int) (locSecond.getStart() - locFirst.getStop());
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
        int result = 0;
        if ( this == that ) {
            result = 0;
        }
        else if(GenomeLoc.isUnmapped(this))
            result = 1;
        else if(GenomeLoc.isUnmapped(that))
            result = -1;
        else {
            final int cmpContig = compareContigs(that);

            if ( cmpContig != 0 ) {
                result = cmpContig;
            } else {
                if ( this.getStart() < that.getStart() ) result = -1;
                if ( this.getStart() > that.getStart() ) result = 1;
            }

            // TODO: and error is being thrown because we are treating reads with the same start positions
            // but different stop as out of order
            //if ( this.getStop() < that.getStop() ) return -1;
            //if ( this.getStop() > that.getStop() ) return 1;
        }

        //System.out.printf("this vs. that = %s %s => %d%n", this, that, result);
        return result;
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
