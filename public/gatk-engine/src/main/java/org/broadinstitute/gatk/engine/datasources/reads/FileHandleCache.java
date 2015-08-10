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

package org.broadinstitute.gatk.engine.datasources.reads;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Caches frequently used  file handles.  Right now, caches only a single file handle.
 * TODO: Generalize to support arbitrary file handle caches.
 */
public class FileHandleCache {
    /**
     * The underlying data structure storing file handles.
     */
    private final FileHandleStorage fileHandleStorage;

    /**
     * How many file handles should be kept open at once.
     */
    private final int cacheSize;

    /**
     * A uniquifier: assign a unique ID to every instance of a file handle.
     */
    private final Map<SAMReaderID,Integer> keyCounter = new HashMap<SAMReaderID,Integer>();

    /**
     * A shared lock, private so that outside users cannot notify it.
     */
    private final Object lock = new Object();

    /**
     * Indicates how many file handles are outstanding at this point.
     */
    private int numOutstandingFileHandles = 0;

    /**
     * Create a new file handle cache of the given cache size.
     * @param cacheSize how many readers to hold open at once.
     */
    public FileHandleCache(final int cacheSize) {
        this.cacheSize = cacheSize;
        fileHandleStorage = new FileHandleStorage();
    }

    /**
     * Retrieves or opens a file handle for the given reader ID.
     * @param key The ke
     * @return A file input stream from the cache, if available, or otherwise newly opened.
     */
    public FileInputStream claimFileInputStream(final SAMReaderID key) {
        synchronized(lock) {
            FileInputStream inputStream = findExistingEntry(key);
            if(inputStream == null) {
                try {
                    // If the cache is maxed out, wait for another file handle to emerge.
                    if(numOutstandingFileHandles >= cacheSize)
                        lock.wait();
                }
                catch(InterruptedException ex) {
                    throw new ReviewedGATKException("Interrupted while waiting for a file handle");
                }
                inputStream = openInputStream(key);
            }
            numOutstandingFileHandles++;

            //System.out.printf("Handing input stream %s to thread %s%n",inputStream,Thread.currentThread().getId());
            return inputStream;
        }
    }

    /**
     * Releases the current reader and returns it to the cache.
     * @param key The reader.
     * @param inputStream The stream being used.
     */
    public void releaseFileInputStream(final SAMReaderID key, final FileInputStream inputStream) {
        synchronized(lock) {
            numOutstandingFileHandles--;
            UniqueKey newID = allocateKey(key);
            fileHandleStorage.put(newID,inputStream);
            // Let any listeners know that another file handle has become available.
            lock.notify();
        }
    }

    /**
     * Finds an existing entry in the storage mechanism.
     * @param key Reader.
     * @return a cached stream, if available.  Otherwise,
     */
    private FileInputStream findExistingEntry(final SAMReaderID key) {
        int existingHandles = getMostRecentUniquifier(key);

        // See if any of the keys currently exist in the repository.
        for(int i = 0; i <= existingHandles; i++) {
            UniqueKey uniqueKey = new UniqueKey(key,i);
            if(fileHandleStorage.containsKey(uniqueKey))
                return fileHandleStorage.remove(uniqueKey);
        }

        return null;
    }

    /**
     * Gets the most recent uniquifier used for the given reader.
     * @param reader Reader for which to determine uniqueness.
     * @return
     */
    private int getMostRecentUniquifier(final SAMReaderID reader) {
        if(keyCounter.containsKey(reader))
            return keyCounter.get(reader);
        else return -1;
    }

    private UniqueKey allocateKey(final SAMReaderID reader) {
        int uniquifier = getMostRecentUniquifier(reader)+1;
        keyCounter.put(reader,uniquifier);
        return new UniqueKey(reader,uniquifier);
    }

    private FileInputStream openInputStream(final SAMReaderID reader) {
        try {
            return new FileInputStream(reader.getSamFilePath());
        }
        catch(IOException ex) {
            throw new GATKException("Unable to open input file");
        }
    }

    private void closeInputStream(final FileInputStream inputStream) {
        try {
            inputStream.close();
        }
        catch(IOException ex) {
            throw new GATKException("Unable to open input file");
        }
    }

    /**
     * Actually contains the file handles, purging them as they get too old.
     */
    private class FileHandleStorage extends LinkedHashMap<UniqueKey,FileInputStream> {
        /**
         * Remove the oldest entry
         * @param entry Entry to consider removing.
         * @return True if the cache size has been exceeded.  False otherwise.
         */
        @Override
        protected boolean removeEldestEntry(Map.Entry<UniqueKey,FileInputStream> entry) {
            synchronized (lock) {
                if(size() > cacheSize) {
                    keyCounter.put(entry.getKey().key,keyCounter.get(entry.getKey().key)-1);
                    closeInputStream(entry.getValue());

                    return true;
                }
            }
            return false;
        }
    }

    /**
     * Uniquifies a key by adding a numerical uniquifier.
     */
    private class UniqueKey {
        /**
         * The file handle's key.
         */
        private final SAMReaderID key;

        /**
         * A uniquifier, so that multiple of the same reader can exist in the cache.
         */
        private final int uniqueID;

        public UniqueKey(final SAMReaderID reader, final int uniqueID) {
            this.key = reader;
            this.uniqueID = uniqueID;
        }

        @Override
        public boolean equals(Object other) {
            if(!(other instanceof UniqueKey))
                return false;
            UniqueKey otherUniqueKey = (UniqueKey)other;
            return key.equals(otherUniqueKey.key) && this.uniqueID == otherUniqueKey.uniqueID;
        }

        @Override
        public int hashCode() {
            return key.hashCode();
        }
    }



}
