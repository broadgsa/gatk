package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;

/**
 * Adapter (decorator) class for rod iterators. The "raw" rod iterator wrapped into this class
 * should be capable of reading the underlying ROD data file and iterating over successive 
 * genomic locations. The purpose of this adapter is to provide additional seekForward() method:
 * upon a call to this method, the decorated iterator will fastforward to the specified position.
 * NOTE 1: if a particular ROD data file is allowed to have multiple records (lines)
 * associated with the same location, the "raw" iterator must be capable of dealing with this situation
 * by loading all such records at once on a call to next().
 * NOTE 2: the object represented by this class is still a unidirectional iterator: after a call to seekForward(),
 * subsequent calls to seekForward() or next() will work from the position the iterator was fastforwarded to.   
 * @author asivache
 *
 * @param <ROD>
 */
public class RODIterator<ROD extends ReferenceOrderedDatum> implements Iterator<ROD> {
    private PushbackIterator<ROD> it;
    private ROD current = null;
    private GenomeLoc position = null;

    public RODIterator(Iterator<ROD> it) {
        this.it = new PushbackIterator<ROD>(it);
    }

    @Override
    public boolean hasNext() { return it.hasNext(); }

    @Override
    public ROD next() {
        ROD next = it.next();
        if( next != null ) {
            position = next.getLocation().clone();
            current = next;
        }
        return next;
    }

//    @Override
//    public boolean hasNext() { return current != null || it.hasNext(); }
//
//    @Override
//    public ROD next() {
//        if ( current != null ) {
//            ROD prev = current;
//            current = null;
//            return prev;
//        } else {
//            ROD next = it.next();
//            if( next != null ) {
//                position = next.getLocation().clone();
//                //current = next;
//            }
//
//            return next;
//        }
//    }

    /**
     * Returns the current position of this iterator.
     * @return Current position of the iterator, or null if no position exists.
     */
    public GenomeLoc position() {
        return position;
    }

    /**
     * Seeks forward in the file until we reach (or cross) a record at contig / pos
     * If we don't find anything and cross beyond contig / pos, we return null;
     * subsequent call to next() will return the first record located after the specified 
     * position in this case. Otherwise, the first ROD record at or overlapping with
     * the specified position is returned; the subsequent call to next() will return the
     * next ROD record. 
     *
     * NOTE 1: the location object <code>loc</code> should be a single point (not an interval);
     * ROD locations, however, can be extended intervals, in which case first ROD that overlaps the specified
     * position will be returned.
     *  
     * NOTE 2: seekForward() is not exactly like next(): if we are strictly past a record, seekForward will not
     * see it, but it will be returning the "current" record (i.e. the record returned by last call to next() or
     * seekForward()) over and over again and will NOT advance the iterator for as long as the current record's location 
     * overlaps with the query position.
     *  
     * @param loc point-like genomic location to fastforward to.
     * @return ROD object at (or overlapping with) the specified position, or null if no such ROD exists.
     */
    public ROD seekForward(final GenomeLoc loc) {
        final boolean DEBUG = false;

        ROD result = null;

        //if (current != null && current.getName().equals("interval")) {
        //    boolean contains = current.getLocation().containsP(loc);
        //    System.out.printf("  %s : current is %s, seeking to %s, contains %b%n", current.getName(), current.getLocation(), loc, contains);
        //}

        if ( current != null && current.getLocation().containsP(loc) )
            return current;                    

        if ( DEBUG ) System.out.printf("  *** starting seek to %s %d (contig %d) from current location %s %d%n", loc.getContig(), loc.getStart(),
        			loc.getContigIndex(),current==null?"null":current.getLocation().getContig(), current==null?-1:current.getLocation().getStart());
        while ( hasNext() ) {
            ROD proposed = next();
            if( proposed == null )
                continue;
            //System.out.printf("    -> Seeking to %s %d AT %s %d%n", contigName, pos, current.getContig(), current.getStart());
            if ( DEBUG ) System.out.println("   proposed at "+proposed.getLocation()+"; contig index="+proposed.getLocation().getContigIndex());
            boolean containedP = proposed.getLocation().containsP(loc);
            //System.out.printf("  %s -> Seeking to %s, at %s => contains = %b%n", current.getName(), loc, current.getLocation(), containedP);
            int cmp = proposed.getLocation().compareTo(loc);
            if ( cmp < 0 ) {
            	if ( DEBUG ) System.out.println("    we are before...");
                // current occurs before loc, continue searching
                continue;
            }
            else if ( cmp == 0 || containedP ) {
            	if ( DEBUG ) System.out.println("    we found overlap...");
                result = proposed;
                break;
            } else {
            	if ( DEBUG ) System.out.println("    we are after...");
                // current is after loc
                it.pushback(proposed);
                break;
            }
        }

        if ( DEBUG ) {
            if ( result != null )
                System.out.printf("    ### Found %s%n", result.getLocation());
        }

        // make a note that the iterator last seeked to the specified position
        current = result;
        position = loc.clone();

        // we ran out of elements or found something
        return result;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}