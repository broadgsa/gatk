package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;

public class RODIterator<ROD extends ReferenceOrderedDatum> implements Iterator<ROD> {
    private PushbackIterator<ROD> it;
    private ROD current = null;
    private GenomeLoc position = null;

    public RODIterator(Iterator<ROD> it) {
        this.it = new PushbackIterator<ROD>(it);
    }

    public boolean hasNext() { return it.hasNext(); }
    public ROD next() {
        ROD next = it.next();
        if( next != null ) {
            position = next.getLocation().clone();
            current = next;
        }
        return next;
    }

    /**
     * Returns the current position of this iterator.
     * @return Current position of the iterator, or null if no position exists.
     */
    public GenomeLoc position() {
        return position;
    }

    /**
     * Seeks forward in the file until we reach (or cross) a record at contig / pos
     * If we don't find anything and cross beyond contig / pos, we return null
     * Otherwise we return the first object who's start is at pos
     *
     * @param loc
     * @return
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