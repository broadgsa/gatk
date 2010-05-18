package org.broadinstitute.sting.utils.file;

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
    private static final String lockString = ".lock";
    private final File lockFile;
    private FileLock lock = null;
    FileChannel fc = null;

    /**
     * create a file system, given a base file to which a lock string gets appended.
     * @param baseLocation the base file location
     */
    public FSLock(File baseLocation) {
        lockFile = new File(baseLocation.getAbsoluteFile() + lockString);
        lockFile.deleteOnExit();
    }

    /**
     * lock the file
     *
     * @return boolean true if we obtained a lock
     */
    public boolean lock() {
        if (lock != null) throw new IllegalStateException("Unable to lock on file " + lockFile + " there is already a lock active");
        if (lockFile.exists()) {
            System.err.println("exits!!"); 
            return false;
        }
        try {
            fc = new RandomAccessFile(lockFile,"rw").getChannel();
            lock = fc.lock();
            return lock != null && lock.isValid();
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to create lock file named " + lockFile,e);
        } catch (IOException e) {
            throw new StingException("Unable to create lock file named " + lockFile,e);
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
