package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 *
 */
public abstract class GenomeLocProcessingTracker {
    /**
     * Information about processing locations and their owners
     */
    public static final class ProcessingLoc implements Comparable<ProcessingLoc> {
        private final GenomeLoc loc;
        private final String owner;

        /**
         * Create a loc that's already owned
         * @param loc
         * @param owner
         */
        public ProcessingLoc(GenomeLoc loc, String owner) {
            if ( loc == null || owner == null ) {
                throw new ReviewedStingException("BUG: invalid ProcessingLoc detected: " + loc + " owner " + owner);
            }

            this.loc = loc;
            this.owner = owner;
        }

        public GenomeLoc getLoc() {
            return loc;
        }

        public String getOwner() {
            return owner;
        }

        public boolean isOwnedBy(String name) {
            return getOwner().equals(name);
        }

        public String toString() { return String.format("ProcessingLoc(%s,%s)", loc, owner); }

        public boolean equals(Object other) {
            if (other instanceof ProcessingLoc )
                return this.loc.equals(((ProcessingLoc)other).loc) && this.owner.equals(((ProcessingLoc)other).owner);
            else
                return false;
        }

        public int compareTo(ProcessingLoc other) {
            return this.getLoc().compareTo(other.getLoc());
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Code to claim intervals for processing and query for their ownership
    //
    // --------------------------------------------------------------------------------
    /**
     * Queries the current database if a location is owned.  Does not guarantee that the
     * loc can be owned in a future call, though.
     *
     * @param loc
     * @return
     */
    public boolean locIsOwned(GenomeLoc loc) {
        return findOwner(loc) != null;
    }

    // in general this isn't true for the list of locs, as they definitely can occur out of order
    protected static ProcessingLoc findOwnerInSortedList(GenomeLoc loc, List<ProcessingLoc> locs) {
        int i = Collections.binarySearch(locs, new ProcessingLoc(loc, "ignore"));
        return i < 0 ? null : locs.get(i);
    }

    protected static ProcessingLoc findOwnerInUnsortedList(GenomeLoc loc, List<ProcessingLoc> locs) {
        for ( ProcessingLoc l : locs ) {
            if ( l.getLoc().equals(loc) )
                return l;
        }

        return null;
    }

    public ProcessingLoc findOwner(GenomeLoc loc) {
        return findOwnerInUnsortedList(loc, getProcessingLocs());
    }

    /**
     * The workhorse routine.  Attempt to claim processing ownership of loc, with my name.
     * This is an atomic operation -- other threads / processes will wait until this function
     * returns.  The return result is the ProcessingLoc object describing who owns this
     * location.  If the location isn't already claimed and we now own the location, the pl owner
     * will be myName.  Otherwise, the name of the owner can found in the pl.
     *
     * @param loc
     * @param myName
     * @return
     */
    public abstract ProcessingLoc claimOwnership(GenomeLoc loc, String myName);

    /**
     * A higher-level, and more efficient, interface to obtain the next location we own.  Takes an
     * iterator producing objects that support the getLocation() interface, and returns the next
     * object in that stream that we can claim ownership of.  Returns null if we run out of elements
     * during the iteration.
     *
     * Can be more efficiently implemented in subclasses to avoid multiple unlocking
     *
     * @param iterator
     * @param myName
     * @return
     */
    public <T extends HasGenomeLocation> T claimOwnershipOfNextAvailable(Iterator<T> iterator, String myName) {
        while ( iterator.hasNext() ) {
            T elt = iterator.next();
            GenomeLoc loc = elt.getLocation();
            ProcessingLoc proc = claimOwnership(loc, myName);

            if ( proc.isOwnedBy(myName) )
                return elt;
            // if not, we continue our search
        }

        // we never found an object, just return it.
        return null;
    }

    /**
     * Returns the list of currently owned locations, updating the database as necessary.
     * DO NOT MODIFY THIS LIST! As with all parallelizing data structures, the list may be
     * out of date immediately after the call returns, or may be updating on the fly.
     *
     * This is really useful for printing, counting, etc. operations that aren't mission critical
     *
     * @return
     */
    protected abstract List<ProcessingLoc> getProcessingLocs();

    protected void close() {
        // by default we don't do anything
    }
}
