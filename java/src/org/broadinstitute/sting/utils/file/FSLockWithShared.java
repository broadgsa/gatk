package org.broadinstitute.sting.utils.file;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.StingException;

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
     * create a file system, given a base file to which a lock string gets appended.
     * @param baseFile File descriptor of file to lock
     */
    public FSLockWithShared(File baseFile) {
        file = baseFile;
    }

    /**
     * Get a shared (read) lock on a file
     * Cannot get shared lock if it does not exist
     * @return boolean true if we obtained a lock
     */
    public boolean sharedLock() {

        // get read-only file channel
        try {
            channel = new RandomAccessFile(file, "r").getChannel();
        }
        catch (IOException e) {
            logger.debug("Unable to lock file (could not open read only file channel)");
            return false;
        }
        // get shared lock (third argument is true)
        try {
            lock = channel.tryLock(0, Long.MAX_VALUE, true);
            if (lock == null) {
                logger.debug("Unable to lock file because there is already a lock active.");
                return false;
            }
        }
        catch (ClosedChannelException e) {
            logger.debug("Unable to lock file because the file channel is closed.");
            return false;
        }
        catch (OverlappingFileLockException e) {
            logger.debug("Unable to lock file because you already have a lock on this file.");
            return false;
        }
        catch (IOException e) {
            logger.debug("Unable to lock file (due to IO exception)");
            return false;
        }
        return true;
    }

    /**
     * Get an exclusive lock on a file
     * @return boolean true if we obtained a lock
     */
    public boolean exclusiveLock() {

        // read/write file channel is necessary for exclusive lock
        try {
            channel = new RandomAccessFile(file, "rw").getChannel();
        }
        catch (Exception e) {
            logger.debug("Unable to lock file (could not open read/write file channel)");
            // do we need to worry about deleting file here? Does RandomAccessFile will only create file if successful?
            return false;
        }

        // get exclusive lock (third argument is false)
        try {
            lock = channel.tryLock(0, Long.MAX_VALUE, false);
            if (lock == null) {
                logger.debug("Unable to lock file because there is already a lock active.");
                return false;
            }
            else return true;
        }
        catch (ClosedChannelException e) {
            logger.debug("Unable to lock file because the file channel is closed.");
            return false;
        }
        catch (OverlappingFileLockException e) {
            logger.debug("Unable to lock file because you already have a lock on this file.");
            return false;
        }
        catch (IOException e) {
            logger.debug("Unable to lock file (due to IO exception)");
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
            throw new StingException("An error occurred while unlocking file", e);
        }
    }
}
