package org.broadinstitute.sting.utils.distributedutils;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;

/**
 * User: depristo
 * Date: 1/19/11
 * Time: 8:24 AM
 *
 * A reentrant lock that supports multi-threaded locking as well as a shared file lock on a common
 * file in the file system.  It itself a shared memory reenterant lock to managed thread safety and
 * contains a SharedFileLock to handle the file integrity.
 */
public class SharedFileThreadSafeLock extends ClosableReentrantLock {
    private static Logger logger = Logger.getLogger(SharedFileThreadSafeLock.class);
    protected static final boolean DEBUG = false;

    private final SharedFileLock fileLock;

    /**
     * Create a SharedFileThreadSafeLock object locking the file
     * @param file
     */
    public SharedFileThreadSafeLock(File file, int nRetries, long milliSecPerTry, int ID) {
        super();
        this.fileLock = new SharedFileLock(file, nRetries, milliSecPerTry, ID);
    }

    public SharedFileThreadSafeLock(File file, int ID) {
        this(file, SharedFileLock.DEFAULT_N_TRIES, SharedFileLock.DEFAULT_MILLISECONDS_PER_TRY, ID);
    }

    @Override
    public void close() {
        super.close();
        fileLock.close();
    }

    @Override
    public int getHoldCount() {
        if ( super.getHoldCount() != fileLock.getHoldCount() )
            throw new ReviewedStingException("BUG: unequal hold counts.  threadlock = " + super.getHoldCount() + ", filelock = " + fileLock.getHoldCount());
        return super.getHoldCount();
    }

    @Override
    public boolean ownsLock() {
        return super.isHeldByCurrentThread() && fileLock.ownsLock();
    }

    /**
     * Two stage [threading then file] locking mechanism.  Reenterant in that multiple lock calls will be
     * unwound appropriately.  Uses file channel lock *after* thread locking.
     */
    @Override
    public void lock() {
        if ( DEBUG ) logger.warn("Attempting SharedFileThreadSafe lock: " + Thread.currentThread().getName());
        if ( DEBUG ) logger.warn("  going for thread lock: " + Thread.currentThread().getName());
        super.lock();
        if ( DEBUG ) logger.warn("  going for file lock: " + Thread.currentThread().getName());
        fileLock.lock(); // todo -- should this be in a try?
    }

    @Override
    public void unlock() {
        if ( DEBUG ) logger.warn("  releasing filelock: " + Thread.currentThread().getName());
        fileLock.unlock();
        if ( DEBUG ) logger.warn("  releasing threadlock: " + Thread.currentThread().getName());
        super.unlock();
        if ( DEBUG ) logger.warn("  unlock() complete: " + Thread.currentThread().getName());
    }
}
