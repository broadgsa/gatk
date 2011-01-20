package org.broadinstitute.sting.utils.threading;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.IOException;
import java.nio.channels.*;
import java.util.concurrent.locks.ReentrantLock;

/**
 * User: depristo
 * Date: 1/19/11
 * Time: 8:24 AM
 *
 * A reentrant lock that supports multi-threaded locking as well as a shared file lock on a common
 * file in the file system.  It itself a shared memory reenterant lock to managed thread safety and a
 * FileChannel FileLock to handle the file integrity
 */
public class SharedFileThreadSafeLock extends ClosableReentrantLock {
    private static Logger logger = Logger.getLogger(SharedFileThreadSafeLock.class);
    private static final boolean DEBUG = false;

    /** The file lock itself that guards the file */
    FileLock fileLock;

    /** the channel object that 'owns' the file lock, and we use to request the lock */
    FileChannel channel;
    int fileLockReentrantCounter = 0;

    /**
     * Create a SharedFileThreadSafeLock object locking the file associated with channel
     * @param channel
     */
    public SharedFileThreadSafeLock(FileChannel channel) {
        super();
        this.channel = channel;
    }

    public void close() {
        try {
            channel.close();
        }
        catch (IOException e) {
            throw new UserException("Count not close channel " + channel, e);
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
            try {
                // Precondition -- lock is always null while we don't have a lock
                if ( fileLock != null )
                    throw new ReviewedStingException("BUG: lock() function called when a lock already is owned!");
                if ( fileLockReentrantCounter == 0 )
                    fileLock = channel.lock();
                fileLockReentrantCounter++;
                if ( DEBUG ) logger.warn("  Have filelock: "  + Thread.currentThread().getName());
            } catch (ClosedChannelException e) {
                throw new ReviewedStingException("Unable to lock file because the file channel is closed. " + channel, e);
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
                fileLock = null;
                if ( DEBUG ) logger.warn("  released filelock: " + Thread.currentThread().getName());
            } else {
                if ( DEBUG ) logger.warn("  skipping filelock release, reenterring unlock via multiple threads " + Thread.currentThread().getName());
            }
        } catch ( IOException e ) {
            throw new ReviewedStingException("Could not free lock on file " + channel, e);
        } finally {
            if ( DEBUG ) logger.warn("  going to release threadlock: " + Thread.currentThread().getName());
            super.unlock();
            if ( DEBUG ) logger.warn("  released threadlock: " + Thread.currentThread().getName());
        }
    }
}
