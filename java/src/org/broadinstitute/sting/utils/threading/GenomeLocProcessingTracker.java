package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.RandomAccessFile;
import java.util.*;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 */
public abstract class GenomeLocProcessingTracker {
    private static Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private Map<GenomeLoc, ProcessingLoc> processingLocs;
    private ClosableReentrantLock lock;
    private long nLockingEvents = 0;

    // --------------------------------------------------------------------------------
    //
    // Factory methods for creating ProcessingTrackers
    //
    // --------------------------------------------------------------------------------

    public static NoOpGenomeLocProcessingTracker createNoOp() {
        return new NoOpGenomeLocProcessingTracker();
    }

    public static SharedMemoryGenomeLocProcessingTracker createSharedMemory() {
        return new SharedMemoryGenomeLocProcessingTracker(new ClosableReentrantLock());
    }

    public static FileBackedGenomeLocProcessingTracker createFileBackedThreaded(File sharedFile, GenomeLocParser parser) {
        return createFileBacked(sharedFile, parser, false);
    }

    public static FileBackedGenomeLocProcessingTracker createFileBackedDistributed(File sharedFile, GenomeLocParser parser) {
        return createFileBacked(sharedFile, parser, true);
    }

    private static FileBackedGenomeLocProcessingTracker createFileBacked(File sharedFile, GenomeLocParser parser, boolean useFileLockToo) {
        try {
            //logger.warn("Creating file backed GLPT at " + sharedFile);
            RandomAccessFile raFile = new RandomAccessFile(sharedFile, "rws");
            ClosableReentrantLock lock = useFileLockToo ? new SharedFileThreadSafeLock(raFile.getChannel()) : new ClosableReentrantLock();
            return new FileBackedGenomeLocProcessingTracker(sharedFile, raFile, parser, lock);
        }
        catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(sharedFile, e);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Creating ProcessingTrackers
    //
    // --------------------------------------------------------------------------------
    public GenomeLocProcessingTracker(ClosableReentrantLock lock) {
        processingLocs = new HashMap<GenomeLoc, ProcessingLoc>();
        this.lock = lock;
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
        // fast path to check if we already have the existing genome loc in memory for ownership claims
        // getProcessingLocs() may be expensive [reading from disk, for example] so we shouldn't call it
        // unless necessary
        ProcessingLoc x = findOwnerInMap(loc, processingLocs);
        return x == null ? findOwnerInMap(loc, updateLocs()) : x;
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
    public ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
        // processingLocs is a shared memory synchronized object, and this
        // method is synchronized, so we can just do our processing
        lock();
        try {
            ProcessingLoc owner = findOwner(loc);

            if ( owner == null ) { // we are unowned
                owner = new ProcessingLoc(loc, myName);
                registerNewLoc(owner);
            }

            return owner;
            //logger.warn(String.format("%s.claimOwnership(%s,%s) => %s", this, loc, myName, owner));
        } finally {
            unlock();
        }
    }

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
        OwnershipIterator<T> myIt = new OwnershipIterator<T>(iterator, myName);
        return myIt.next();
    }

    public <T extends HasGenomeLocation> Iterable<T> onlyOwned(Iterator<T> iterator, String myName) {
        return new OwnershipIterator<T>(iterator, myName);
    }

    protected class OwnershipIterator<T extends HasGenomeLocation> implements Iterator<T>, Iterable<T> {
        Iterator<T> subit;
        String myName;

        public OwnershipIterator(Iterator<T> subit, String myName) {
            this.subit = subit;
            this.myName = myName;
        }

        /**
         * Will return true for all elements of subit, even if we can't get ownership of some of the future
         * elements and so will return null there
         * @return
         */
        public boolean hasNext() {
            return subit.hasNext();
        }

        /**
         * High performance iterator that only locks and unlocks once per claimed object found.  Avoids
         * locking / unlocking for each query
         *
         * @return an object of type T owned by this thread, or null if none of the remaining object could be claimed
         */
        public T next() {
            lock();
            try {
                while ( subit.hasNext() ) {
                    T elt = subit.next();
                    //logger.warn("Checking elt for ownership " + elt);
                    GenomeLoc loc = elt.getLocation();
                    ProcessingLoc proc = claimOwnership(loc, myName);

                    if ( proc.isOwnedBy(myName) )
                        return elt;
                    // if not, we continue our search
                }

                // we never found an object, just return it.
                return null;
            } finally {
                unlock();
            }
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        public Iterator<T> iterator() {
            return this;
        }
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
    protected Collection<ProcessingLoc> getProcessingLocs() {
        return updateLocs().values();
    }

    private Map<GenomeLoc, ProcessingLoc> updateLocs() {
        lock();
        try {
            for ( ProcessingLoc p : readNewLocs() )
                processingLocs.put(p.getLocation(), p);
            return processingLocs;
        } finally {
            unlock();
        }
    }


    // --------------------------------------------------------------------------------
    //
    // Low-level accessors / manipulators and utility functions
    //
    // --------------------------------------------------------------------------------

    private final void lock() {
        if ( ! lock.isHeldByCurrentThread() )
            nLockingEvents++;
        lock.lock();
    }

    private final void unlock() {
        lock.unlock();
    }

    protected static ProcessingLoc findOwnerInCollection(GenomeLoc loc, Collection<ProcessingLoc> locs) {
        for ( ProcessingLoc l : locs ) {
            if ( l.getLocation().equals(loc) )
                return l;
        }

        return null;
    }

    protected static ProcessingLoc findOwnerInMap(GenomeLoc loc, Map<GenomeLoc,ProcessingLoc> locs) {
        return locs.get(loc);
    }


    // --------------------------------------------------------------------------------
    //
    // Code to override to change the dynamics of the the GenomeLocProcessingTracker
    //
    // --------------------------------------------------------------------------------

    protected void close() {
        lock.close();
        logger.warn("Locking events: " + nLockingEvents);
        // by default we don't do anything
    }

    protected abstract void registerNewLoc(ProcessingLoc loc);
    protected abstract Collection<ProcessingLoc> readNewLocs();
}
