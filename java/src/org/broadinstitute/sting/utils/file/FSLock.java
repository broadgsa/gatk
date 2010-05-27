package org.broadinstitute.sting.utils.file;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;

/**
 * a quick implementation of a file based lock, using the Java NIO classes
 */
public class FSLock {
    // connect to the logger
    private static Logger logger = Logger.getLogger(FSLock.class);

    // the lock string
    private static final String lockString = ".lock";

    // the file we're attempting to lock
    private final File lockFile;

    // the file lock
    private FileLock lock = null;

    // the file channel we open
    FileChannel fc = null;

    /**
     * create a file system, given a base file to which a lock string gets appended.
     * @param baseLocation the base file location
     */
    public FSLock(File baseLocation) {
        lockFile = new File(baseLocation.getAbsoluteFile() + lockString);
        if (lockFile != null) lockFile.deleteOnExit();
    }

    /**
     * lock the file
     *
     * @return boolean true if we obtained a lock
     */
    public synchronized boolean lock() {
        if (lock != null) throw new IllegalStateException("Unable to lock on file " + lockFile + " there is already a lock active");

        // if the file already exists, or we can't write to the parent directory, return false
        if (lockFile.exists() || !lockFile.getParentFile().canWrite())
            return false;
        try {
            fc = new RandomAccessFile(lockFile,"rw").getChannel();
            lock = fc.tryLock();
            return lock != null && lock.isValid();
        } catch (FileNotFoundException e) {
            logger.info("Unable to create lock file (due to no file found exception), file = " + lockFile,e);
            return false;
        } catch (IOException e) {
            logger.info("Unable to create lock file (due to IO exception), file = " + lockFile,e);
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
            if (fc != null)
                fc.close();
            if (lockFile.exists())
                lockFile.delete();
        } catch (IOException e) {
            throw new StingException("Unable to create lock file named " + lockFile,e);
        }
    }
}
