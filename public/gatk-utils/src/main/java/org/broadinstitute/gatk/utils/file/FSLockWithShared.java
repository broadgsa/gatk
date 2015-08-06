/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.file;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.*;
import java.util.concurrent.*;

/**
 * a quick implementation of a file based lock, using the Java NIO classes
 */
public class FSLockWithShared {
    // connect to the logger
    private final static Logger logger = Logger.getLogger(FSLockWithShared.class);

    // the file we're attempting to lock
    private final File file;

    // the file lock
    private FileLock lock = null;

    // the file channel we open
    private FileChannel channel = null;

    // Timeout (in milliseconds) before we give up during non-blocking lock-acquisition calls.
    // Necessary because these "non-blocking" calls can hang if there's a problem with the
    // OS file locking support.
    private int lockAcquisitionTimeout;

    // Default value for lockAcquisitionTimeout when none is explicitly provided
    public static final int DEFAULT_LOCK_ACQUISITION_TIMEOUT_IN_MILLISECONDS = 30 * 1000;

    // Amount of time to wait when trying to shut down the lock-acquisition thread before giving up
    public static final int THREAD_TERMINATION_TIMEOUT_IN_MILLISECONDS = 30 * 1000;

    /**
     * Create a lock associated with the specified File. Use the default lock
     * acquisition timeout of 30 seconds.
     *
     * @param file file to lock
     */
    public FSLockWithShared( final File file ) {
        this.file = file;
        lockAcquisitionTimeout = DEFAULT_LOCK_ACQUISITION_TIMEOUT_IN_MILLISECONDS;
    }

    /**
     * Create a lock associated with the specified File, and set a custom lock
     * acquisition timeout.
     *
     * @param file file to lock
     * @param lockAcquisitionTimeout maximum number of milliseconds to wait during non-blocking
     *                               lock acquisition calls before concluding that there's a
     *                               problem with the OS file locking support and throwing an error.
     */
    public FSLockWithShared( final File file, final int lockAcquisitionTimeout ) {
        this.file = file;
        this.lockAcquisitionTimeout = lockAcquisitionTimeout;
    }

    /**
     * Get a shared (read) lock on a file. Does not block, and returns immediately
     * under normal conditions with the result of the lock acquisition attempt. Will
     * throw an exception if there's a problem with the OS file locking support.
     *
     * @return boolean true if we obtained a lock, false if we failed to obtain one
     */
    public boolean sharedLock() {
        return acquireLockWithTimeout(true);
    }

    /**
     * Get an exclusive (read-write) lock on a file. Does not block, and returns immediately
     * under normal conditions with the result of the lock acquisition attempt. Will
     * throw an exception if there's a problem with the OS file locking support.
     *
     * @return boolean true if we obtained a lock, false if we failed to obtain one
     */
    public boolean exclusiveLock() {
        return acquireLockWithTimeout(false);
    }

    /**
     * Attempt to acquire a lock of the specified type on the file in a background thread.
     * Uses non-blocking lock-acquisition calls that should return immediately, but may
     * get stuck if there's a problem with the OS file locking support. If the call gets
     * stuck and the timeout elapses, throws a UserException, since it's not safe to
     * proceed with a stuck lock acquisition thread (and there's no way to reliably
     * interrupt it once the underlying system call hangs).
     *
     * @param acquireSharedLock if true, request a shared lock rather than an exclusive lock
     * @return true if a lock was acquired, false if we failed
     */
    private boolean acquireLockWithTimeout( final boolean acquireSharedLock ) {
        // Use daemon threads so that hopelessly stuck lock acquisition threads won't prevent the JVM from exiting
        final ExecutorService executor = Executors.newSingleThreadExecutor(new ThreadFactory() {
                                                                               public Thread newThread( Runnable r ) {
                                                                                   Thread lockAcquisitionThread = new Thread(r);
                                                                                   lockAcquisitionThread.setDaemon(true);
                                                                                   return lockAcquisitionThread;
                                                                               }
                                                                           });
        final FutureTask<Boolean> lockAcquisitionTask = new FutureTask<Boolean>(new LockAcquisitionTask(acquireSharedLock));
        boolean lockAcquired = false;

        try {
            executor.execute(lockAcquisitionTask);

            // Wait at most lockAcquisitionTimeout milliseconds for the lock acquisition task to finish.
            lockAcquired = lockAcquisitionTask.get(lockAcquisitionTimeout, TimeUnit.MILLISECONDS);
        }
        // Lock acquisition timeout elapsed. Since we're using NON-BLOCKING lock-acquisition calls,
        // this implies that there's a problem with the OS locking daemon, or locks are not supported.
        // Since it's not safe to proceed with a potentially stuck lock acquisition thread, we need to
        // shut down the JVM in order to kill it.
        catch ( TimeoutException e ) {
            throw new UserException.FileSystemInabilityToLockException(
                    String.format("Timeout of %d milliseconds was reached while trying to acquire a lock on file %s. " +
                                  "Since the GATK uses non-blocking lock acquisition calls that are not supposed to wait, " +
                                  "this implies a problem with the file locking support in your operating system.",
                                  lockAcquisitionTimeout, file.getAbsolutePath()));
        }
        // Lock acquisition thread threw an exception. Need to unpack it via e.getCause()
        catch ( ExecutionException e ) {
            logger.warn(String.format("WARNING: Unable to lock file %s because exception %s occurred with error message %s",
                                      file.getAbsolutePath(),
                                      e.getCause() != null ? e.getCause().getClass().getSimpleName() : "unknown",
                                      e.getCause() != null ? e.getCause().getMessage() : "none"));
            lockAcquired = false;
        }
        // Interrupted while waiting for the lock acquisition thread -- not likely to happen
        catch ( InterruptedException e ) {
            logger.warn(String.format("WARNING: interrupted while attempting to acquire a lock for file %s", file.getAbsolutePath()));
            lockAcquired = false;
        }
        catch ( Exception e ) {
            logger.warn(String.format("WARNING: error while attempting to acquire a lock for file %s. Error message: %s",
                                      file.getAbsolutePath(), e.getMessage()));
            lockAcquired = false;
        }

        shutdownLockAcquisitionTask(executor);

        // Upon failure to acquire a lock, we always call unlock() to close the FileChannel if it was opened
        // and to deal with very hypothetical edge cases where a lock might actually have been acquired despite the
        // lock acquisition thread returning false.
        if ( ! lockAcquired ) {
            unlock();
        }

        return lockAcquired;
    }

    /**
     * Ensures that the lock acquisition task running in the provided executor has cleanly terminated.
     * Throws a UserException if unable to shut it down within the period defined by the THREAD_TERMINATION_TIMEOUT.
     *
     * @param executor ExecutorService executing the lock-acquisition thread
     */
    private void shutdownLockAcquisitionTask( final ExecutorService executor ) {
        boolean shutdownAttemptSucceeded;

        try {
            executor.shutdownNow();
            shutdownAttemptSucceeded = executor.awaitTermination(THREAD_TERMINATION_TIMEOUT_IN_MILLISECONDS, TimeUnit.MILLISECONDS);
        }
        catch ( InterruptedException e ) {
            shutdownAttemptSucceeded = false;
        }

        if ( ! shutdownAttemptSucceeded ) {
            throw new UserException(String.format("Failed to terminate lock acquisition thread while trying to lock file %s. " +
                                                  "Exiting because it's not safe to proceed with this run of the GATK.",
                                                  file.getAbsolutePath()));
        }
    }

    /**
     * Background task that attempts to acquire a lock of the specified type, and returns a boolean
     * indicating success/failure. Uses a non-blocking tryLock() call that should return immediately
     * (but may get stuck if there's a problem with the OS locking daemon).
     */
    private class LockAcquisitionTask implements Callable<Boolean> {
        private final boolean acquireSharedLock;

        public LockAcquisitionTask( final boolean acquireSharedLock ) {
            this.acquireSharedLock = acquireSharedLock;
        }

        public Boolean call() {
            // Get a read-only or read-write file channel, depending on the type of lock
            try {
                channel = new RandomAccessFile(file, acquireSharedLock ? "r" : "rw").getChannel();
            }
            catch ( IOException e ) {
                logger.warn(String.format("WARNING: Unable to lock file %s because we could not open a file channel", file.getAbsolutePath()));
                return false;
            }

            boolean lockAcquired = false;

            try {
                // Non-blocking lock-acquisition call, should return right away. If it doesn't return immediately
                // due to problems with the OS locking daemon, it will potentially be timed-out and interrupted.
                lock = channel.tryLock(0, Long.MAX_VALUE, acquireSharedLock);
                lockAcquired = lock != null;
            }
            catch ( AsynchronousCloseException e ) {
                logger.warn(String.format("WARNING: Unable to lock file %s because the file channel was closed by another thread", file.getAbsolutePath()));
                lockAcquired = false;
            }
            catch ( ClosedChannelException e ) {
                logger.warn(String.format("WARNING: Unable to lock file %s because the file channel is closed.", file.getAbsolutePath()));
                lockAcquired = false;
            }
            catch ( OverlappingFileLockException e ) {
                logger.warn(String.format("WARNING: Unable to lock file %s because you already have a lock on this file.", file.getAbsolutePath()));
                lockAcquired = false;
            }
            catch ( FileLockInterruptionException e ) {
                logger.warn(String.format("WARNING: Interrupted while attempting to lock file %s", file.getAbsolutePath()));
                lockAcquired = false;
            }
            catch ( IOException e ) {
                logger.warn(String.format("WARNING: Unable to lock file %s because an IOException occurred with message: %s.", file.getAbsolutePath(), e.getMessage()));
                lockAcquired = false;
            }

            return lockAcquired;
        }
    }

    /**
     * Unlock the file
     *
     * note: this allows unlocking a file that failed to lock (no required user checks on null locks).
     */
    public void unlock() {
        releaseLock();
        closeChannel();
    }

    private void releaseLock() {
        try {
            if ( lock != null )
                lock.release();
        }
        catch ( ClosedChannelException e ) {
            // if the channel was already closed we don't have to worry
        }
        catch ( IOException e ) {
            throw new UserException(String.format("An error occurred while releasing the lock for file %s", file.getAbsolutePath()), e);
        }
    }

    private void closeChannel() {
        try {
            if ( channel != null )
                channel.close();
        }
        catch ( IOException e ) {
            throw new UserException(String.format("An error occurred while closing channel for file %s", file.getAbsolutePath()), e);
        }
    }
}
