package org.broadinstitute.sting.utils.distributedutils;

import org.apache.log4j.Logger;
import org.apache.lucene.store.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;

/**
 * User: depristo
 * Date: 1/19/11
 * Time: 8:24 AM
 *
 * A reentrant lock for a shared file common file in the file system.  Relies on a a Lucene SimpleFSLock
 * to manage on disk file locking.
 */
public class SharedFileLock extends ClosableReentrantLock { // todo -- kinda gross inheritance.  The super lock is never used
    private static Logger logger = Logger.getLogger(SharedFileLock.class);

    private static final String VERIFY_HOST = System.getProperty("verify.host", "gsa1");
    private static final boolean VERIFY = false;
    private static final int VERIFY_PORT = 5050;

    // 5 minutes => 360 seconds of trying -> failure
    protected static final int DEFAULT_N_TRIES = 1000;
    protected static final long DEFAULT_MILLISECONDS_PER_TRY = 360;

    /** The file we are locking */
    private final File file;

    private final LockFactory lockFactory;
    private Lock fileLock = null;

    /**
     * A counter that indicates the number of 'locks' on this file.
     * If locks == 2, then two unlocks are required
     * before any resources are freed.
     */
    int fileLockReentrantCounter = 0;

    // type of locking
    private final int nRetries;
    private final long milliSecPerTry;

    /**
     * Create a SharedFileThreadSafeLock object locking the file
     * @param file
     */
    public SharedFileLock(File file, int nRetries, long milliSecPerTry, int ID) {
        super();
        this.file = file;
        this.nRetries = nRetries;
        this.milliSecPerTry = milliSecPerTry;

        File lockDir = new File(file.getParent() == null ? "./" : file.getParent());
        try {
            LockFactory factory = new SimpleFSLockFactory(lockDir);
            if ( VERIFY ) { // don't forget to start up the VerifyLockServer
                this.lockFactory = new VerifyingLockFactory((byte)ID, factory, VERIFY_HOST, VERIFY_PORT);
            } else {
                this.lockFactory = factory;
            }
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(lockDir, "Could not create coordination file locking directory " + lockDir, e);
        }
    }

    public SharedFileLock(File file, int ID) {
        this(file, DEFAULT_N_TRIES, DEFAULT_MILLISECONDS_PER_TRY, ID);
    }

    @Override
    public void close() {
        if ( ownsLock() ) throw new ReviewedStingException("closing SharedFileLock while still owned: ownership count " + fileLockReentrantCounter);
    }

    @Override
    public int getHoldCount() {
        return fileLockReentrantCounter;
    }

    @Override
    public boolean ownsLock() {
        return fileLockReentrantCounter > 0;
    }

    // ------------------------------------------------------------------------------------------
    //
    // workhorse routines -- acquiring file locks
    //
    // ------------------------------------------------------------------------------------------

    private boolean obtainFileLock() throws IOException {
        // annoying bug work around for verifylockserver
        if ( VERIFY )
            try {
                return fileLock.obtain(1);
            } catch ( LockObtainFailedException e ) {
                return false;
            }
        else
            return fileLock.obtain();
    }

    /**
     * Two stage [threading then file] locking mechanism.  Reenterant in that multiple lock calls will be
     * unwound appropriately.  Uses file channel lock *after* thread locking.
     */
    @Override
    public void lock() {
        if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            lock() " + Thread.currentThread().getName() + ", fileLockReentrantCounter = " + fileLockReentrantCounter);
        if ( fileLockReentrantCounter++ == 0 ) {
            // Precondition -- lock is always null while we don't have a lock
            if ( fileLock != null )
                throw new ReviewedStingException("BUG: lock() function called when a lock already is owned!");

            int i = 1;
            fileLock = lockFactory.makeLock(file.getName() + ".lock");
            try {
                boolean obtained = obtainFileLock(); // todo -- maybe use intrinsic lock features
                for ( ; ! obtained && i < nRetries; i++ ) {
                    try {
                        //logger.warn("tryLock failed on try " + i + ", waiting " + milliSecPerTry + " millseconds for retry");
                        Thread.sleep(milliSecPerTry);
                    } catch ( InterruptedException e ) {
                        throw new UserException("SharedFileThreadSafeLock interrupted during wait for file lock", e);
                    }
                    obtained = obtainFileLock(); // gross workaround for error in verify server
                }

                if ( i > 1 ) logger.warn("tryLock required " + i + " tries before completing, waited " + i * milliSecPerTry + " millseconds");

                if ( ! obtained ) {
                    fileLock = null;
                    // filelock == null -> we never managed to acquire the lock!
                    throw new UserException("SharedFileThreadSafeLock failed to obtain the lock after " + nRetries + " attempts");
                }

                if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            lock() " + Thread.currentThread().getName() + ", obtained = " + obtained + ", tries = " + i);
            } catch (IOException e) {
                fileLock = null;
                throw new ReviewedStingException("Coordination file could not be created because a lock could not be obtained.", e);
            }
        }
    }

    @Override
    public void unlock() {
        // update for reentrant unlocking
        if ( fileLock == null ) throw new ReviewedStingException("BUG: file lock is null -- file lock was not obtained");
        if ( fileLockReentrantCounter <= 0 ) throw new ReviewedStingException("BUG: file lock counter < 0");

        // this unlock counts as 1 unlock.  If this is our last unlock, actually do something
        if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", count = " + fileLockReentrantCounter);
        if ( --fileLockReentrantCounter == 0 ) {
            try {
                if ( ! fileLock.isLocked() ) throw new ReviewedStingException("BUG: call to unlock() when we don't have a valid lock!");
                fileLock.release();
                if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", actually releasing");
            } catch ( IOException e ) {
                throw new ReviewedStingException("Could not free file lock on file " + file, e);
            } finally {  // make sure we null out the filelock, regardless of our state
                fileLock = null;
            }
        } else {
            if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", skipping, count = " + fileLockReentrantCounter);
        }
    }
}
