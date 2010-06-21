/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.FastaSequenceIndex;
import org.broadinstitute.sting.utils.fasta.FastaSequenceIndexBuilder;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import net.sf.picard.sam.CreateSequenceDictionary;
import org.broadinstitute.sting.utils.file.FSLock;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;

/**
 * Loads reference data from fasta file
 * Looks for fai and dict files, and tries to create them if they don't exist
 */
public class ReferenceDataSource implements ReferenceDataSourceProgressListener {
    private IndexedFastaSequenceFile index;

    /** our log, which we want to capture anything from this class */
    protected static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(ReferenceDataSource.class);

    /**
     * Create reference data source from fasta file
     * @param fastaFile Fasta file to be used as reference
     */
    public ReferenceDataSource(File fastaFile) {
        File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");
        File dictFile = new File(fastaFile.getAbsolutePath().replace(".fasta", ".dict"));

        /*
         if index file does not exist, create it manually
          */
        if (!indexFile.exists()) {
            logger.info(String.format("Index file %s does not exist. Trying to create it now.", indexFile.getAbsolutePath()));
            FileChannel indexChannel;
            FileLock indexLock;
            try {
                // get exclusive lock
                indexChannel = new RandomAccessFile(indexFile, "rw").getChannel();
                if ((indexLock = indexChannel.tryLock(0, Long.MAX_VALUE, false)) == null)
                    throw new StingException("Index file could not be written because a lock could not be obtained." +
                            "If you are running multiple instances of GATK, another process is probably creating this " +
                            "file now. Please wait until it is finished and try again.");
                FastaSequenceIndexBuilder faiBuilder = new FastaSequenceIndexBuilder(fastaFile, this);
                FastaSequenceIndex sequenceIndex = faiBuilder.createIndex();
                FastaSequenceIndexBuilder.saveAsFaiFile(sequenceIndex, indexFile);
                // unlock
                try {
                    if (indexLock != null)
                        indexLock.release();
                    if (indexChannel != null)
                        indexChannel.close();
                }
                catch (Exception e) {
                    throw new StingException("An error occurred while unlocking file:" + indexFile.getAbsolutePath(), e);
                }
            }
            catch (Exception e) {
                throw new StingException("Index file does not exist and could not be created. See error below.", e);
            }
        }

        /*
        * If dict file doesn't exist, try to create it using Picard's CreateSequenceDictionary
        * Currently, dictionary cannot be created without running CreateSequenceDictionary's main routine, hence the
        * argument string
        * This has been filed in trac as (PIC-370) Want programmatic interface to CreateSequenceDictionary
        */
        if (!dictFile.exists()) {
            logger.info(String.format("Index file %s does not exist. Trying to create it now.", indexFile.getAbsolutePath()));
            FileChannel dictChannel;
            FileLock dictLock;
            try {
                // get exclusive lock
                dictChannel = new RandomAccessFile(indexFile, "rw").getChannel();
                if ((dictLock = dictChannel.tryLock(0, Long.MAX_VALUE, false)) == null)

                    throw new StingException("Dictionary file could not be written because a lock could not be obtained." +
                            "If you are running multiple instances of GATK, another process is probably creating this " +
                            "file now. Please wait until it is finished and try again.");

                // create dictionary by calling main routine. Temporary fix - see comment above.
                String args[] = {String.format("r=%s", fastaFile.getAbsolutePath()),
                        String.format("o=%s", dictFile.getAbsolutePath())};
                new CreateSequenceDictionary().instanceMain(args);
                // unlock
                try {
                    if (dictLock != null)
                        dictLock.release();
                    if (dictChannel != null)
                        dictChannel.close();
                }
                catch (Exception e) {
                    throw new StingException("An error occurred while unlocking file:" + indexFile.getAbsolutePath(), e);
                }
            }
            catch (Exception e) {
                throw new StingException("Dictionary file does not exist and could not be created. See error below.", e);
            }
        }

        /*
         * Read reference data by creating an IndexedFastaSequenceFile.
         * A note about thread safety: IndexFastaSequenceFile reads the fasta using dictionary and index files. It will
         * fail if either does not exist, but not if either is currently being written (in which case it exists
         * but is incomplete). To avoid this, obtain shared locks on both files before creating IndexedFastaSequenceFile.
         */

        FileChannel dictChannel;
        FileChannel indexChannel;
        FileLock dictLock;
        FileLock indexLock;
        try {
            // set up dictionary and index locks
            // channel is read only and lock is shared (third argument is true) 
            dictChannel = new RandomAccessFile(dictFile, "r").getChannel();
            if ((dictLock = dictChannel.tryLock(0, Long.MAX_VALUE, true)) == null) {
                throw new StingException("Could not open dictionary file because a lock could not be obtained.");
            }
            indexChannel = new RandomAccessFile(indexFile, "r").getChannel();
            if ((indexLock = indexChannel.tryLock(0, Long.MAX_VALUE, true)) == null) {
                throw new StingException("Could not open dictionary file because a lock could not be obtained.");
            }

            index = new IndexedFastaSequenceFile(fastaFile);

            // unlock/close
            try {
                if (dictLock != null)
                    dictLock.release();
                if (dictChannel != null)
                    dictChannel.close();
            }
            catch (Exception e) {
                throw new StingException("An error occurred while unlocking file:" + dictFile.getAbsolutePath(), e);
            }
            try {
                if (indexLock != null)
                    indexLock.release();
                if (indexChannel != null)
                    indexChannel.close();
            }
            catch (Exception e) {
                throw new StingException("An error occurred while unlocking file:" + indexFile.getAbsolutePath(), e);
            }
        }
        catch (Exception e) {
            throw new StingException(String.format("Error reading fasta file %s. See stack trace below.", fastaFile.getAbsolutePath()), e);
        }
    }

    /**
     * Get indexed fasta file
     * @return IndexedFastaSequenceFile that was created from file
     */
    public IndexedFastaSequenceFile getReference() {
        return this.index;
    }

    /**
     * Notify user of progress in creating fai file
     * @param percent Percent of fasta file read as a percent
     */
    public void percentProgress(int percent) {
        System.out.println(String.format("PROGRESS UPDATE: file is %d percent complete", percent));
    }

}
