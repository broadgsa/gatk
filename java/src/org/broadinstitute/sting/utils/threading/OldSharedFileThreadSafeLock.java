package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.*;

/**
 * User: depristo
 * Date: 1/19/11
 * Time: 8:24 AM
 *
 * A reentrant lock that supports multi-threaded locking as well as a shared file lock on a common
 * file in the file system.  It itself a shared memory reenterant lock to managed thread safety and a
 * FileChannel FileLock to handle the file integrity
 */
public class OldSharedFileThreadSafeLock extends ClosableReentrantLock {
    private static Logger logger = Logger.getLogger(OldSharedFileThreadSafeLock.class);
    private static final boolean DEBUG = false;

    // 100 seconds of trying -> failure
    private static final int DEFAULT_N_TRIES = 1000;
    private static final long DEFAULT_MILLISECONDS_PER_TRY = 100;

    /** The file we are locking */
    private final File file;

    /** The file lock itself that guards the file */
    FileLock fileLock;

    /** the channel object that 'owns' the file lock, and we use to request the lock */
    FileChannel channel;

    /**
     * A counter that indicates the number of 'locks' on this file.
     * If locks == 2, then two unlocks are required
     * before any resources are freed.
     */
    int fileLockReentrantCounter = 0;

    // type of locking
    private final boolean blockOnLock;
    private final int nRetries;
    private final long milliSecPerTry;

    /**
     * Create a SharedFileThreadSafeLock object locking the file
     * @param file
     */
    public OldSharedFileThreadSafeLock(File file, boolean blockOnLock, int nRetries, long milliSecPerTry) {
        super();
        this.file = file;
        this.blockOnLock = blockOnLock;
        this.nRetries = nRetries;
        this.milliSecPerTry = milliSecPerTry;
    }

    public OldSharedFileThreadSafeLock(File file, boolean blockOnLock) {
        this(file, blockOnLock, DEFAULT_N_TRIES, DEFAULT_MILLISECONDS_PER_TRY);
    }


    private FileChannel getChannel() {
        if ( DEBUG ) logger.warn("    Get channel: " + Thread.currentThread().getName() + " channel = " + channel);
        if ( channel == null ) {
            try {
                if ( DEBUG ) logger.warn("    opening channel: " + Thread.currentThread().getName());
                this.channel = new RandomAccessFile(file, "rw").getChannel();
                if ( DEBUG ) logger.warn("    opened channel: " + Thread.currentThread().getName());
            } catch (FileNotFoundException e) {
                throw new UserException.CouldNotCreateOutputFile(file, e);
            }
        }

        return this.channel;
    }

    private void closeChannel() {
        try {
            if ( channel != null ) {
                channel.close();
                channel = null;
            }
        }
        catch (IOException e) {
            throw new UserException("Count not close channel associated with file" + file, e);
        }
    }

    public void close() {
        super.close();
        closeChannel();
    }

    public boolean ownsLock() {
        return super.isHeldByCurrentThread() && fileLockReentrantCounter > 0;
    }

    // ------------------------------------------------------------------------------------------
    //
    // workhorse routines -- acquiring file locks
    //
    // ------------------------------------------------------------------------------------------

    private void acquireFileLock() {
        try {
            // Precondition -- lock is always null while we don't have a lock
            if ( fileLock != null )
                throw new ReviewedStingException("BUG: lock() function called when a lock already is owned!");

            if ( blockOnLock ) {
                //
                // blocking code
                //
                fileLock = getChannel().lock();
            } else {
                //
                // polling code
                //
                int i = 0;
                for ( ; fileLock == null && i < nRetries; i++ ) {
                    fileLock = getChannel().tryLock();
                    if ( fileLock == null ) {
                        try {
                            //logger.warn("tryLock failed on try " + i + ", waiting " + milliSecPerTry + " millseconds for retry");
                            Thread.sleep(milliSecPerTry);
                        } catch ( InterruptedException e ) {
                            throw new UserException("SharedFileThreadSafeLock interrupted during wait for file lock", e);
                        }
                    }
                }
                if ( i > 1 ) logger.warn("tryLock required " + i + " tries before completing, waited " + i * milliSecPerTry + " millseconds");

                if ( fileLock == null ) {
                    // filelock == null -> we never managed to acquire the lock!
                    throw new UserException("SharedFileThreadSafeLock failed to obtain the lock after " + nRetries + " attempts");
                }
            }
            if ( DEBUG ) logger.warn("  Have filelock: "  + Thread.currentThread().getName());
        } catch (ClosedChannelException e) {
            throw new ReviewedStingException("Unable to lock file because the file channel is closed. " + file, e);
        } catch (FileLockInterruptionException e) {
            throw new ReviewedStingException("File lock interrupted", e);
        } catch (NonWritableChannelException e) {
            throw new ReviewedStingException("File channel not writable", e);
        } catch (OverlappingFileLockException e) {
            // this only happens when multiple threads are running, and one is waiting
            // for the lock above and we come here.
            throw new ReviewedStingException("BUG: Failed to acquire lock, should never happen.");
        } catch (IOException e) {
            throw new ReviewedStingException("Coordination file could not be created because a lock could not be obtained.", e);
        }
    }

    /**
     * Two stage [threading then file] locking mechanism.  Reenterant in that multiple lock calls will be
     * unwound appropriately.  Uses file channel lock *after* thread locking.
     */
    @Override
    public void lock() {
        if ( DEBUG ) logger.warn("Attempting threadlock: " + Thread.currentThread().getName());

        if ( super.isHeldByCurrentThread() ) {
            if ( DEBUG ) logger.warn("  Already have threadlock, continuing: " + Thread.currentThread().getName());
            super.lock();                          // call the lock here so we can call unlock later
            fileLockReentrantCounter++;            // inc. the file lock counter
            return;
        } else {
            super.lock();
            if ( DEBUG ) logger.warn("  Have thread-lock, going for filelock: " + Thread.currentThread().getName());
            if ( fileLockReentrantCounter == 0 )
                acquireFileLock();
            fileLockReentrantCounter++;
        }
    }

    @Override
    public void unlock() {
        try {
            // update for reentrant unlocking
            fileLockReentrantCounter--;
            if ( fileLockReentrantCounter < 0 ) throw new ReviewedStingException("BUG: file lock counter < 0");

            if ( fileLock != null && fileLockReentrantCounter == 0 ) {
                if ( ! fileLock.isValid() ) throw new ReviewedStingException("BUG: call to unlock() when we don't have a valid lock!");

                if ( DEBUG ) logger.warn("  going to release filelock: " + Thread.currentThread().getName());
                fileLock.release();
                closeChannel();
                fileLock = null;
                if ( DEBUG ) logger.warn("  released filelock: " + Thread.currentThread().getName());
            } else {
                if ( DEBUG ) logger.warn("  skipping filelock release, reenterring unlock via multiple threads " + Thread.currentThread().getName());
            }
        } catch ( IOException e ) {
            throw new ReviewedStingException("Could not free lock on file " + file, e);
        } finally {
            if ( DEBUG ) logger.warn("  going to release threadlock: " + Thread.currentThread().getName());
            super.unlock();
            if ( DEBUG ) logger.warn("  released threadlock: " + Thread.currentThread().getName());
        }
    }
}
