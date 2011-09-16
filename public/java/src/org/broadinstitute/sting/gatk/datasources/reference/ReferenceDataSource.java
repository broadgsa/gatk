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

package org.broadinstitute.sting.gatk.datasources.reference;

import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.FastaSequenceIndexBuilder;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.sam.CreateSequenceDictionary;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.broadinstitute.sting.utils.file.FileSystemInabilityToLockException;

import java.io.File;

/**
 * Loads reference data from fasta file
 * Looks for fai and dict files, and tries to create them if they don't exist
 */
public class ReferenceDataSource {
    private IndexedFastaSequenceFile index;

    /** our log, which we want to capture anything from this class */
    protected static org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(ReferenceDataSource.class);

    /**
     * Create reference data source from fasta file
     * @param fastaFile Fasta file to be used as reference
     */
    public ReferenceDataSource(File fastaFile) {

        // does the fasta file exist? check that first...
        if (!fastaFile.exists())
            throw new UserException("The fasta file you specified (" + fastaFile.getAbsolutePath() + ") does not exist.");

        File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");
        File dictFile;
        if (fastaFile.getAbsolutePath().endsWith("fa")) {
            dictFile = new File(fastaFile.getAbsolutePath().replace(".fa", ".dict"));
        }
        else
         dictFile = new File(fastaFile.getAbsolutePath().replace(".fasta", ".dict"));

        /*
         if index file does not exist, create it manually
          */
        if (!indexFile.exists()) {
            logger.info(String.format("Index file %s does not exist. Trying to create it now.", indexFile.getAbsolutePath()));
            FSLockWithShared indexLock = new FSLockWithShared(indexFile,true);
            try {
                // get exclusive lock
                if (!indexLock.exclusiveLock())
                    throw new UserException.CouldNotCreateReferenceIndexFileBecauseOfLock(dictFile);
                FastaSequenceIndexBuilder faiBuilder = new FastaSequenceIndexBuilder(fastaFile, true);
                FastaSequenceIndex sequenceIndex = faiBuilder.createIndex();
                FastaSequenceIndexBuilder.saveAsFaiFile(sequenceIndex, indexFile);
            }
            catch(FileSystemInabilityToLockException ex) {
                logger.info("Unable to create write lock: " + ex.getMessage());
                logger.info("Skipping index creation.");
            }
            catch(UserException e) {
                // Rethrow all user exceptions as-is; there should be more details in the UserException itself. 
                throw e;    
            }
            catch (Exception e) {
                // If lock creation succeeded, the failure must have been generating the index.
                // If lock creation failed, just skip over index creation entirely.
                throw new UserException.CouldNotCreateReferenceIndexFile(indexFile, e);
            }
            finally {
                indexLock.unlock();
            }
        }

        /*
        * If dict file doesn't exist, try to create it using Picard's CreateSequenceDictionary
        * Currently, dictionary cannot be created without running CreateSequenceDictionary's main routine, hence the
        * argument string
        * This has been filed in trac as (PIC-370) Want programmatic interface to CreateSequenceDictionary
        */
        if (!dictFile.exists()) {
            logger.info(String.format("Dict file %s does not exist. Trying to create it now.", dictFile.getAbsolutePath()));

            /*
             * Please note another hack here: we have to create a temporary file b/c CreateSequenceDictionary cannot
             * create a dictionary file if that file is locked.
             */

            // get read lock on dict file so nobody else can read it
            FSLockWithShared dictLock = new FSLockWithShared(dictFile,true);
            try {
                // get shared lock on dict file so nobody else can start creating it
                if (!dictLock.exclusiveLock())
                    throw new UserException.CouldNotCreateReferenceIndexFileBecauseOfLock(dictFile);
                // dict will be written to random temporary file in same directory (see note above)
                File tempFile = File.createTempFile("dict", null, dictFile.getParentFile());
                tempFile.deleteOnExit();

                // create dictionary by calling main routine. Temporary fix - see comment above.
                String args[] = {String.format("r=%s", fastaFile.getAbsolutePath()),
                        String.format("o=%s", tempFile.getAbsolutePath())};
                new CreateSequenceDictionary().instanceMain(args);

                if (!tempFile.renameTo(dictFile))
                    throw new UserException("Error transferring temp file " + tempFile + " to dict file " + dictFile);
            }
            catch(FileSystemInabilityToLockException ex) {
                logger.info("Unable to create write lock: " + ex.getMessage());
                logger.info("Skipping dictionary creation.");
            }
            catch (Exception e) {
                // If lock creation succeeded, the failure must have been generating the index.
                // If lock creation failed, just skip over index creation entirely.
                throw new UserException.CouldNotCreateReferenceIndexFile(dictFile, e);
            }
            finally {
                dictLock.unlock();
            }
        }

        /*
         * Read reference data by creating an IndexedFastaSequenceFile.
         * A note about thread safety: IndexFastaSequenceFile reads the fasta using dictionary and index files. It will
         * fail if either does not exist, but not if either is currently being written (in which case it exists
         * but is incomplete). To avoid this, obtain shared locks on both files before creating IndexedFastaSequenceFile.
         */

        FSLockWithShared dictLock = new FSLockWithShared(dictFile,true);
        FSLockWithShared indexLock = new FSLockWithShared(indexFile,true);
        try {
            try {
                if (!dictLock.sharedLock()) {
                    throw new ReviewedStingException("Could not open dictionary file because a lock could not be obtained.");
                }
            }
            catch(FileSystemInabilityToLockException ex) {
                logger.info(String.format("Unable to create a lock on dictionary file: %s",ex.getMessage()));
                logger.info("Treating existing dictionary file as complete.");
            }

            try {
                if (!indexLock.sharedLock()) {
                    throw new ReviewedStingException("Could not open index file because a lock could not be obtained.");
                }
            }
            catch(FileSystemInabilityToLockException ex) {
                logger.info(String.format("Unable to create a lock on index file: %s",ex.getMessage()));
                logger.info("Treating existing index file as complete.");
            }

            index = new CachingIndexedFastaSequenceFile(fastaFile);

        } catch (IllegalArgumentException e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, "Could not read reference sequence.  The FASTA must have either a .fasta or .fa extension", e);
        }
        catch (Exception e) {
            throw new UserException.CouldNotReadInputFile(fastaFile, e);
        }
        finally {
            dictLock.unlock();
            indexLock.unlock();
        }
    }

    /**
     * Get indexed fasta file
     * @return IndexedFastaSequenceFile that was created from file
     */
    public IndexedFastaSequenceFile getReference() {
        return this.index;
    }
}
