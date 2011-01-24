package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.locks.ReentrantLock;

/**
 *
 */
public abstract class GenomeLocProcessingTracker {
    private final static Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private final static SimpleDateFormat STATUS_FORMAT = new SimpleDateFormat("HH:mm:ss,SSS");
    private final static int DEFAULT_OWNERSHIP_ITERATOR_SIZE = 100;

    private final static String GOING_FOR_LOCK = "going_for_lock";
    private final static String HAVE_LOCK = "have_lock";
    private final static String RUNNING = "running";

    private final Map<GenomeLoc, ProcessingLoc> processingLocs;
    private final ClosableReentrantLock lock;
    private final PrintStream status;

    protected SimpleTimer writeTimer = new SimpleTimer("writeTimer");
    protected SimpleTimer readTimer = new SimpleTimer("readTimer");
    protected SimpleTimer lockWaitTimer = new SimpleTimer("lockWaitTimer");
    protected long nLocks = 0, nWrites = 0, nReads = 0;

    // --------------------------------------------------------------------------------
    //
    // Factory methods for creating ProcessingTrackers
    //
    // --------------------------------------------------------------------------------

    public static GenomeLocProcessingTracker createNoOp() {
        return new NoOpGenomeLocProcessingTracker();
    }

    public static GenomeLocProcessingTracker createSharedMemory() {
        return new SharedMemoryGenomeLocProcessingTracker(new ClosableReentrantLock());
    }

    public static GenomeLocProcessingTracker createFileBackedThreaded(File sharedFile, GenomeLocParser parser, PrintStream status) {
        return createFileBacked(sharedFile, parser, false, false, status);
    }

    public static GenomeLocProcessingTracker createFileBackedDistributed(File sharedFile, GenomeLocParser parser, boolean blockingP, PrintStream status) {
        return createFileBacked(sharedFile, parser, blockingP, true, status);
    }

    private static FileBackedGenomeLocProcessingTracker createFileBacked(File sharedFile, GenomeLocParser parser, boolean blockP, boolean useFileLockToo, PrintStream status) {
        //logger.warn("Creating file backed GLPT at " + sharedFile);
        ClosableReentrantLock lock = useFileLockToo ? new SharedFileThreadSafeLock(sharedFile, blockP) : new ClosableReentrantLock();
        return new FileBackedGenomeLocProcessingTracker(sharedFile, parser, lock, status);
    }

    // --------------------------------------------------------------------------------
    //
    // Creating ProcessingTrackers
    //
    // --------------------------------------------------------------------------------
    public GenomeLocProcessingTracker(ClosableReentrantLock lock, PrintStream status) {
        this.processingLocs = new HashMap<GenomeLoc, ProcessingLoc>();
        this.status = status;
        this.lock = lock;
        printStatusHeader();
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
    public final boolean locIsOwned(GenomeLoc loc, String id) {
        return findOwner(loc, id) != null;
    }

    protected final ProcessingLoc findOwner(GenomeLoc loc, String id) {
        // fast path to check if we already have the existing genome loc in memory for ownership claims
        // getProcessingLocs() may be expensive [reading from disk, for example] so we shouldn't call it
        // unless necessary
        ProcessingLoc x = findOwnerInMap(loc, processingLocs);
        return x == null ? findOwnerInMap(loc, updateLocs(id)) : x;
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
    public final ProcessingLoc claimOwnership(GenomeLoc loc, String myName) {
        // processingLocs is a shared memory synchronized object, and this
        // method is synchronized, so we can just do our processing
        lock(myName);
        try {
            ProcessingLoc owner = findOwner(loc, myName);

            if ( owner == null ) { // we are unowned
                owner = new ProcessingLoc(loc, myName);
                registerNewLocsWithTimers(Arrays.asList(owner), myName);
            }

            return owner;
            //logger.warn(String.format("%s.claimOwnership(%s,%s) => %s", this, loc, myName, owner));
        } finally {
            unlock(myName);
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
    public final <T extends HasGenomeLocation> T claimOwnershipOfNextAvailable(Iterator<T> iterator, String myName) {
        OwnershipIterator<T> myIt = new OwnershipIterator<T>(iterator, myName, 1);
        return myIt.next();
    }

    public final <T extends HasGenomeLocation> Iterable<T> onlyOwned(Iterator<T> iterator, String myName) {
        return new OwnershipIterator<T>(iterator, myName);
    }

    protected final class OwnershipIterator<T extends HasGenomeLocation> implements Iterator<T>, Iterable<T> {
        private final Iterator<T> subit;
        private final String myName;
        private final Queue<T> cache;
        private final int cacheSize;

        public OwnershipIterator(Iterator<T> subit, String myName) {
            this(subit, myName, DEFAULT_OWNERSHIP_ITERATOR_SIZE);
        }

        public OwnershipIterator(Iterator<T> subit, String myName, int cacheSize) {
            this.subit = subit;
            this.myName = myName;
            cache = new LinkedList<T>();
            this.cacheSize = cacheSize;
        }

        /**
         * Will return true for all elements of subit, even if we can't get ownership of some of the future
         * elements and so will return null there
         * @return
         */
        public final boolean hasNext() {
            return cache.peek() != null || subit.hasNext();
        }

        /**
         * High performance iterator that only locks and unlocks once per claimed object found.  Avoids
         * locking / unlocking for each query
         *
         * @return an object of type T owned by this thread, or null if none of the remaining object could be claimed
         */
        public final T next() {
            T elt = cache.poll();
            if ( elt != null)
                return elt;
            else {
                // cache is empty, we need to fill up the cache and return the first element of the queue
                lock(myName);
                try {

                    // read once the database of owners at the start
                    updateLocs(myName);

                    boolean done = false;
                    Queue<ProcessingLoc> pwns = new LinkedList<ProcessingLoc>(); // ;-)
                    while ( !done && cache.size() < cacheSize && subit.hasNext() ) {
                        elt = subit.next();
                        //logger.warn("Checking elt for ownership " + elt);
                        GenomeLoc loc = elt.getLocation();

                        ProcessingLoc owner = findOwnerInMap(loc, processingLocs);

                        if ( owner == null ) { // we are unowned
                            owner = new ProcessingLoc(loc, myName);
                            pwns.offer(owner);
                            if ( ! cache.offer(elt) ) throw new ReviewedStingException("Cache offer unexpectedly failed");
                            if ( GenomeLoc.isUnmapped(loc) ) done = true;
                        }
                        // if not, we continue our search
                    }

                    registerNewLocsWithTimers(pwns, myName);

                    // we've either filled up the cache or run out of elements.  Either way we return
                    // the first element of the cache. If the cache is empty, we return null here.
                    //logger.warn("Cache size is " + cache.size());
                    //logger.warn("Cache contains " + cache);

                    return cache.poll();
                } finally {
                    unlock(myName);
                }
            }
        }

        public final void remove() {
            throw new UnsupportedOperationException();
        }

        public final Iterator<T> iterator() {
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
    protected final Collection<ProcessingLoc> getProcessingLocs(String myName) {
        return updateLocs(myName).values();
    }

    private final Map<GenomeLoc, ProcessingLoc> updateLocs(String myName) {
        lock(myName);
        try {
            readTimer.restart();
            for ( ProcessingLoc p : readNewLocs() )
                processingLocs.put(p.getLocation(), p);
            readTimer.stop();
            nReads++;
            return processingLocs;
        } finally {
            unlock(myName);
        }
    }

    protected final void registerNewLocsWithTimers(Collection<ProcessingLoc> plocs, String myName) {
        writeTimer.restart();
        registerNewLocs(plocs);
        nWrites++;
        writeTimer.stop();
    }

    // --------------------------------------------------------------------------------
    //
    // Low-level accessors / manipulators and utility functions
    //
    // --------------------------------------------------------------------------------
    private final boolean hasStatus() {
        return status != null;
    }

    private final void printStatusHeader() {
        if ( hasStatus() ) status.printf("process.id\thr.time\ttime\tstate%n");
    }

    private final void printStatus(String id, long machineTime, String state) {
        // prints a line like processID human-readable-time machine-time state
        if ( hasStatus() ) {
            status.printf("%s\t%s\t%d\t%s%n", id, STATUS_FORMAT.format(machineTime), machineTime, state);
            status.flush();
        }
    }

    private final void lock(String id) {
        lockWaitTimer.restart();
        boolean hadLock = lock.ownsLock();
        if ( ! hadLock ) {
            nLocks++;
            printStatus(id, lockWaitTimer.currentTime(), GOING_FOR_LOCK);
        }
        lock.lock();
        lockWaitTimer.stop();
        if ( ! hadLock ) printStatus(id, lockWaitTimer.currentTime(), HAVE_LOCK);
    }

    private final void unlock(String id) {
        lock.unlock();
        if ( ! lock.ownsLock() ) printStatus(id, lockWaitTimer.currentTime(), RUNNING);
    }

    protected final static ProcessingLoc findOwnerInCollection(GenomeLoc loc, Collection<ProcessingLoc> locs) {
        for ( ProcessingLoc l : locs ) {
            if ( l.getLocation().equals(loc) )
                return l;
        }

        return null;
    }

    protected final static ProcessingLoc findOwnerInMap(GenomeLoc loc, Map<GenomeLoc,ProcessingLoc> locs) {
        return locs.get(loc);
    }

    // useful code for getting
    public final long getNLocks() { return nLocks; }
    public final long getNReads() { return nReads; }
    public final long getNWrites() { return nWrites; }
    public final double getTimePerLock() { return lockWaitTimer.getElapsedTime() / Math.max(nLocks, 1); }
    public final double getTimePerRead() { return readTimer.getElapsedTime() / Math.max(nReads,1); }
    public final double getTimePerWrite() { return writeTimer.getElapsedTime() / Math.max(nWrites,1); }

    // --------------------------------------------------------------------------------
    //
    // Code to override to change the dynamics of the the GenomeLocProcessingTracker
    //
    // --------------------------------------------------------------------------------

    protected void close() {
        lock.close();
        if ( hasStatus() ) status.close();
        //logger.warn("Locking events: " + nLocks);
    }

    protected abstract void registerNewLocs(Collection<ProcessingLoc> plocs);
    protected abstract Collection<ProcessingLoc> readNewLocs();
}
