package org.broadinstitute.sting.utils.file;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.ClosedChannelException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;

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

    /**
     * A bit of experimental code for Siva at Partners.  Conditionally throw an
     * exception in the case where an unknown failure occurs, in an effort to stave
     * off disabled nfs file locks.
     */
    private boolean throwExceptionOnUnknownFailure = false;

    /**
     * create a file system, given a base file to which a lock string gets appended.
     * @param baseFile File descriptor of file to lock
     */
    public FSLockWithShared(File baseFile) {
        file = baseFile;
    }

    public FSLockWithShared(File baseFile,boolean throwExceptionOnUnknownFailure) {
        this(baseFile);
        this.throwExceptionOnUnknownFailure = throwExceptionOnUnknownFailure;
    }

    /**
     * Get a shared (read) lock on a file
     * Cannot get shared lock if it does not exist
     * @return boolean true if we obtained a lock
     * @throws FileSystemInabilityToLockException in cases of unexpected failure to capture lock.
     */
    public boolean sharedLock() throws FileSystemInabilityToLockException {

        // get read-only file channel
        try {
            channel = new RandomAccessFile(file, "r").getChannel();
        }
        catch (IOException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s (could not open read only file channel)",file.getAbsolutePath()));
            return false;
        }
        // get shared lock (third argument is true)
        try {
            lock = channel.tryLock(0, Long.MAX_VALUE, true);
            if (lock == null) {
                logger.warn(String.format("WARNING: Unable to lock file %s because there is already a lock active.",file.getAbsolutePath()));
                return false;
            }
        }
        catch (ClosedChannelException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s because the file channel is closed.",file.getAbsolutePath()));
            return false;
        }
        catch (OverlappingFileLockException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s because you already have a lock on this file.",file.getAbsolutePath()));
            return false;
        }
        catch (IOException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s: %s.",file.getAbsolutePath(),e.getMessage()));
            if(throwExceptionOnUnknownFailure)
                throw new FileSystemInabilityToLockException(e.getMessage(),e);
            else
                return false;
        }
        return true;
    }

    /**
     * Get an exclusive lock on a file
     * @return boolean true if we obtained a lock
     * @throws FileSystemInabilityToLockException in cases of unexpected failure to capture lock.
     */
    public boolean exclusiveLock() throws FileSystemInabilityToLockException {

        // read/write file channel is necessary for exclusive lock
        try {
            channel = new RandomAccessFile(file, "rw").getChannel();
        }
        catch (Exception e) {
            logger.warn(String.format("WARNING: Unable to lock file %s (could not open read/write file channel)",file.getAbsolutePath()));
            // do we need to worry about deleting file here? Does RandomAccessFile will only create file if successful?
            return false;
        }

        // get exclusive lock (third argument is false)
        try {
            lock = channel.tryLock(0, Long.MAX_VALUE, false);
            if (lock == null) {
                logger.warn(String.format("WARNING: Unable to lock file %s because there is already a lock active.",file.getAbsolutePath()));
                return false;
            }
            else return true;
        }
        catch (ClosedChannelException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s because the file channel is closed.",file.getAbsolutePath()));
            return false;
        }
        catch (OverlappingFileLockException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s because you already have a lock on this file.",file.getAbsolutePath()));
            return false;
        }
        catch (IOException e) {
            logger.warn(String.format("WARNING: Unable to lock file %s: %s.",file.getAbsolutePath(),e.getMessage()));
            if(throwExceptionOnUnknownFailure)
                throw new FileSystemInabilityToLockException(e.getMessage(),e);
            else
                return false;
        }
    }
    
    /**
     * unlock the file
     *
     * note: this allows unlocking a file that failed to lock (no required user checks on null locks).
     */
    public void unlock() {
        try {
            if (lock != null)
                lock.release();
            if (channel != null)
                channel.close();
        }
        catch (Exception e) {
            throw new ReviewedStingException("An error occurred while unlocking file", e);
        }
    }
}
