package org.broadinstitute.sting.utils.threading;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 1/13/11
 * Time: 9:38 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class GenomeLocProcessingTracker {
    /**
     * Information about processing locations and their owners
     */
    public static final class ProcessingLoc {
        GenomeLoc loc;
        String owner;

        /**
         * Create an unowned loc
         * @param loc
         */
        public ProcessingLoc(GenomeLoc loc) {
            this(loc, null);

        }

        /**
         * Create a loc that's already owned
         * @param loc
         * @param owner
         */
        public ProcessingLoc(GenomeLoc loc, String owner) {
            this.loc = loc;
            this.owner = owner;
        }

        public GenomeLoc getLoc() {
            return loc;
        }

        public String getOwner() {
            return owner;
        }

        public void setOwner(String owner) {
            this.owner = owner;
        }

        public boolean isOwned() { return owner != null; }
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

    public ProcessingLoc findOwner(GenomeLoc loc) {
        for ( ProcessingLoc l : getProcessingLocs() ) {
            if ( l.getLoc().equals(loc) )
                return l;
        }

        return null;
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
     * Returns the list of currently owned locations, updating the database as necessary.
     * DO NOT MODIFY THIS LIST! As with all parallelizing data structures, the list may be
     * out of date immediately after the call returns, or may be updating on the fly.
     *
     * This is really useful for printing, counting, etc. operations that aren't mission critical
     *
     * @return
     */
    protected abstract List<ProcessingLoc> getProcessingLocs();
}
