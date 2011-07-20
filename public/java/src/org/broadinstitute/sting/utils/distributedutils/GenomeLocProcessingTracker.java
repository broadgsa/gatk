package org.broadinstitute.sting.utils.distributedutils;

import net.sf.picard.reference.IndexedFastaSequenceFile;
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
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Abstract base class to coordinating data processing by a collecting for processes / threads.
 *
 * Conceptually, the genome is viewed as a collection of non-overlapping genome location:
 *
 * chr1:1-10
 * chr1:11-20
 * chr1:21-30
 * etc.
 *
 * This class, and it's concrete derived classes, provide the ability to claim individual locations
 * as "mine", and exclude other processes / threads from processing them.  At the lowest-level this
 * is implemented by the claimOwnership(loc, name) function, that returns true if loc free (unclaimed)
 * and makes name the owner of loc.  High-level, and more efficient operations provide claiming
 * iterators over streams of objects implementing the HasGenomeLocation interface, so that you can
 * write code that looks like:
 *
 * for ( GenomeLoc ownedLoc : onlyOwned(allLocsToProcess.iterator) ) {
 *   doSomeWork(ownedLoc)
 *
 * Much of the code in this class is actually surrounding debugging and performance metrics code.
 * The actual synchronization code is separated out into the ClosableReentrantLock() system
 * and the two abstract functions:
 *
 *  protected abstract void registerNewLocs(Collection<ProcessingLoc> plocs);
 *  protected abstract Collection<ProcessingLoc> readNewLocs();
 *
 * That maintain the state of the tracker.
 *
 * That is, the ProcessingTracker is made of two components: a thread / process locking system and
 * a subclass that implements the methods to record new claimed state changes and to read out updates
 * that may have occurred by another thread or process.
 *
 * NOTE: this class assumes that all threads / processes are working with the same set of potential
 * GenomeLocs to own.  Claiming chr1:1-10 and then chr1:5-6 is allowed by the system.  Basically,
 * you only can stake claim to GenomeLocs that are .equal().
 */
public abstract class GenomeLocProcessingTracker {
    private final static Logger logger = Logger.getLogger(FileBackedGenomeLocProcessingTracker.class);
    private final static SimpleDateFormat STATUS_FORMAT = new SimpleDateFormat("HH:mm:ss,SSS");
    private final static int DEFAULT_OWNERSHIP_ITERATOR_SIZE = 1;

    /**
     * Useful state strings for printing status
     */
    private final static String GOING_FOR_LOCK = "going_for_lock";
    private final static String RELEASING_LOCK = "releasing_lock";
    private final static String HAVE_LOCK = "have_lock";
    private final static String RUNNING = "running";

    /**
     * A map, for efficiency, that allows quick lookup of the processing loc for a
     * given GenomeLoc.  The map points from loc -> loc / owner as a ProcessingLoc
     */
    private final Map<GenomeLoc, ProcessingLoc> processingLocs;

    /**
     * The locking object used to protect data from simulatanous access by multiple
     * threads or processes.
     */
    private final ClosableReentrantLock lock;

    /** A stream for writing status messages.  Can be null if we aren't writing status */
    private final PrintStream status;

    //
    // Timers for recording performance information
    // Note -- these cannot be used because this class isn't thread safe, and neither are the
    // timers, so they result in invalid operations w.r.t. the SimpleTimer contract
    //
//    protected final SimpleTimer writeTimer = new SimpleTimer("writeTimer");
//    protected final SimpleTimer readTimer = new SimpleTimer("readTimer");
//    protected final SimpleTimer lockWaitTimer = new SimpleTimer("lockWaitTimer");
    protected final SimpleTimer timer = new SimpleTimer();
    protected long nLocks = 0, nWrites = 0, nReads = 0;

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
    // Code to override to change the dynamics of the the GenomeLocProcessingTracker
    //
    // --------------------------------------------------------------------------------

    protected void close() {
        lock.close();
        if ( status != null ) status.close();
    }

    /**
     * Takes a collection of newly claimed (i.e., previous unclaimed) genome locs
     * and the name of their owner and "registers" this data in some persistent way that's
     * visible to all threads / processes communicating via this GenomeLocProcessingTracker.
     *
     * Could be a in-memory data structure (a list) if we are restricting ourselves to intra-memory
     * parallelism, a locked file on a shared file system, or a server we communicate with.
     *
     * @param plocs
     */
    protected abstract void registerNewLocs(Collection<ProcessingLoc> plocs);

    /**
     * The inverse of the registerNewLocs() function.  Looks at the persistent data store
     * shared by all threads / processes and returns the ones that have appeared since the last
     * call to readNewLocs().  Note that we expect the pair of registerNewLocs and readNewLocs to
     * include everything, even locs registered by this thread / process.  For example:
     *
     * readNewLocs() => List()
     * registerNewLocs(List(x, y,)) => void
     * readNewLocs() => List(x,y))
     *
     * even for this thread or process.
     * @return
     */
    protected abstract Collection<ProcessingLoc> readNewLocs();


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
    public final ProcessingLoc claimOwnership(final GenomeLoc loc, final String myName) {
        // processingLocs is a shared memory synchronized object, and this
        // method is synchronized, so we can just do our processing
        return new WithLock<ProcessingLoc>(myName) {
            public ProcessingLoc doBody() {
                ProcessingLoc owner = findOwner(loc, myName);
                if ( owner == null ) { // we are unowned
                    owner = new ProcessingLoc(loc, myName);
                    registerNewLocsWithTimers(Arrays.asList(owner), myName);
                }
                return owner;
            }
        }.run();
    }


    // --------------------------------------------------------------------------------
    //
    // High-level iterator-style interface to claiming ownership
    //
    // --------------------------------------------------------------------------------

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

    private final class OwnershipIterator<T extends HasGenomeLocation> implements Iterator<T>, Iterable<T> {
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
            if ( cache.peek() != null)
                return cache.poll();
            else {
                // cache is empty, we need to fill up the cache and return the first element of the queue
                return new WithLock<T>(myName) {
                    public T doBody() {
                        // read once the database of owners at the start
                        updateAndGetProcessingLocs(myName);

                        boolean done = false;
                        Queue<ProcessingLoc> pwns = new LinkedList<ProcessingLoc>(); // ;-)
                        while ( !done && cache.size() < cacheSize && subit.hasNext() ) {
                            final T elt = subit.next();
                            GenomeLoc loc = elt.getLocation();

                            ProcessingLoc owner = processingLocs.get(loc);

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
                        return cache.poll();
                    }
                }.run();
            }
        }

        public final void remove() {
            throw new UnsupportedOperationException();
        }

        public final Iterator<T> iterator() {
            return this;
        }
    }

    // --------------------------------------------------------------------------------
    //
    // private / protected low-level accessors / manipulators and utility functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Useful debugging function that returns the ProcessingLoc who owns loc.  ID
     * is provided for debugging purposes
     * @param loc
     * @param id
     * @return
     */
    protected final ProcessingLoc findOwner(GenomeLoc loc, String id) {
        // fast path to check if we already have the existing genome loc in memory for ownership claims
        // getProcessingLocs() may be expensive [reading from disk, for example] so we shouldn't call it
        // unless necessary
        ProcessingLoc x = processingLocs.get(loc);
        return x == null ? updateAndGetProcessingLocs(id).get(loc) : x;
    }

    /**
     * Returns the list of currently owned locations, updating the database as necessary.
     * DO NOT MODIFY THIS MAP! As with all parallelizing data structures, the list may be
     * out of date immediately after the call returns, or may be updating on the fly.
     * @return
     */
    protected final Map<GenomeLoc, ProcessingLoc> updateAndGetProcessingLocs(String myName) {
        return new WithLock<Map<GenomeLoc, ProcessingLoc>>(myName) {
            public Map<GenomeLoc, ProcessingLoc> doBody() {
//                readTimer.restart();
                for ( ProcessingLoc p : readNewLocs() )
                    processingLocs.put(p.getLocation(), p);
//                readTimer.stop();
                nReads++;
                return processingLocs;
            }
        }.run();
    }

    /**
     * Wrapper around registerNewLocs that also times the operation
     *
     * @param plocs
     * @param myName
     */
    protected final void registerNewLocsWithTimers(Collection<ProcessingLoc> plocs, String myName) {
//        writeTimer.restart();
        registerNewLocs(plocs);
        nWrites++;
//        writeTimer.stop();
    }

    private final void printStatusHeader() {
        if ( status != null ) status.printf("process.id\thr.time\ttime\tstate%n");
    }

    private final void printStatus(String id, long machineTime, String state) {
        // prints a line like processID human-readable-time machine-time state
        if ( status != null  ) {
            status.printf("%s\t%s\t%d\t%s%n", id, STATUS_FORMAT.format(machineTime), machineTime, state);
            status.flush();
        }
    }


    /**
     * Lock the data structure, preventing other threads / processes from reading and writing to the
     * common store
     * @param id the name of the process doing the locking
     */
    private final void lock(String id) {
        //lockWaitTimer.restart();
        boolean hadLock = lock.ownsLock();
        if ( ! hadLock ) {
            nLocks++;
            //printStatus(id, lockWaitTimer.currentTime(), GOING_FOR_LOCK);
        }
        lock.lock();
        //lockWaitTimer.stop();
        //if ( ! hadLock ) printStatus(id, lockWaitTimer.currentTime(), HAVE_LOCK);
    }

    /**
     * Unlock the data structure, allowing other threads / processes to read and write to the common store
     * @param id the name of the process doing the unlocking
     */
    private final void unlock(String id) {
        if ( lock.getHoldCount() == 1 ) printStatus(id, timer.currentTime(), RELEASING_LOCK);
        lock.unlock();
        if ( ! lock.ownsLock() ) printStatus(id, timer.currentTime(), RUNNING);
    }

    // useful code for getting
    public final long getNLocks() { return nLocks; }
    public final long getNReads() { return nReads; }
    public final long getNWrites() { return nWrites; }
//    public final double getTimePerLock() { return lockWaitTimer.getElapsedTime() / Math.max(nLocks, 1); }
//    public final double getTimePerRead() { return readTimer.getElapsedTime() / Math.max(nReads,1); }
//    public final double getTimePerWrite() { return writeTimer.getElapsedTime() / Math.max(nWrites,1); }

    // --------------------------------------------------------------------------------
    //
    // Java-style functional form for with lock do { x };
    //
    // --------------------------------------------------------------------------------

    /**
     * Private utility class that executes doBody() method with the lock() acquired and
     * handles property unlock()ing the system, even if an error occurs.  Allows one to write
     * clean code like:
     *
     * new WithLock<Integer>(name) {
     *   public Integer doBody() { doSomething(); return 1; }
     * }.run()
     *
     * @param <T> the return type of the doBody() method
     */
    private abstract class WithLock<T> {
        private final String myName;

        public WithLock(String myName) {
            this.myName = myName;
        }

        protected abstract T doBody();

        public T run() {
            boolean locked = false;
            try {
                lock(myName);
                locked = true;
                return doBody();
            } finally {
                if (locked) unlock(myName);
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // main function for testing performance
    //
    // --------------------------------------------------------------------------------
    public static void main(String[] args) {
        //BasicConfigurator.configure();

        final String ref = args[0];
        final File file = new File(args[1]);
        final int cycles = Integer.valueOf(args[2]);

        File referenceFile = new File(ref);
        try {
            final IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(referenceFile);
            final String chr1 = fasta.getSequenceDictionary().getSequence(1).getSequenceName();
            final GenomeLocParser genomeLocParser = new GenomeLocParser(fasta);

            final class MyTest {
                String name;
                GenomeLocProcessingTracker tracker;

                MyTest(String name, GenomeLocProcessingTracker tracker) {
                    this.name = name;
                    this.tracker = tracker;
                }

                public void execute(int cycles) {
                    SimpleTimer delta = new SimpleTimer("delta");
                    SimpleTimer timer = new SimpleTimer("none");

                    if ( file.exists() ) file.delete();
                    timer.start();
                    delta.start();
                    for ( int i = 1; i < cycles; i++ ) {
                        tracker.claimOwnership(genomeLocParser.createGenomeLoc(chr1, i, i+1), "ABCDEFGHIJKL");
                        if ( i % 1000 == 0 ) {
                            System.out.printf("%s\t%d\t%d\t%.4f\t%.4f%n", name, i, timer.currentTime(), timer.getElapsedTime(), delta.getElapsedTime() );
                            delta.restart();
                        }
                    }
                }
            }

            System.out.printf("name\tcycle\tcurrent.time\telapsed.time\tdelta%n");
            new MyTest("in-memory", new SharedMemoryGenomeLocProcessingTracker(new ClosableReentrantLock())).execute(cycles);
            new MyTest("nio", new FileBackedGenomeLocProcessingTracker(file, genomeLocParser, new ClosableReentrantLock(), null)).execute(cycles);
            new MyTest("nio-file-lock", new FileBackedGenomeLocProcessingTracker(file, genomeLocParser, new SharedFileThreadSafeLock(file,1), null)).execute(cycles);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }
}
