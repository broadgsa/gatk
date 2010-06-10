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

import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.FastaSequenceIndex;
import org.broadinstitute.sting.utils.fasta.FastaSequenceIndexBuilder;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import net.sf.picard.sam.CreateSequenceDictionary;

import java.io.File;
import java.io.IOException;
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
     *
     * Please note: This should be re-done when PIC-370 is fixed. (See below.) At that time, we may consider loading
     * sequenceIndex and sequenceDict here, instead of in IndexedFastaSequenceFile. 
     */
    public ReferenceDataSource(File fastaFile) {
        File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");
        File dictFile = new File(fastaFile.getAbsolutePath().replace(".fasta", ".dict"));
        FastaSequenceIndex sequenceIndex; // stores FastaSequenceIndex if file doesn't exist

        /*
        Note: this code is temporary. See commented code below, which will be used when thread safety is resolved.
         */
        if (!indexFile.exists()) {
            FastaSequenceIndexBuilder faiBuilder = new FastaSequenceIndexBuilder(fastaFile, this);
            sequenceIndex = faiBuilder.sequenceIndex;
            try {
                index = new IndexedFastaSequenceFile(fastaFile, sequenceIndex);
            }
            catch (Exception e) {
                throw new StingException("Could not load fasta file. Stack trace below. ", e);
            }
        }
        else {
            try {
                index = new IndexedFastaSequenceFile(fastaFile);
            }
            catch (Exception e) {
                throw new StingException("Could not load fasta file. Stack trace below.", e);
            }
        }

        /*
        * If dict file doesn't exist, try to create it using Picard's CreateSequenceDictionary
        * Currently, dictionary cannot be created without running CreateSequenceDictionary's main routine, hence the
        * argument string
        * This has been filed in trac as (PIC-370) Want programmatic interface to CreateSequenceDictionary
        *
        * Todo: Currently working on making the creation of dictionary and index thread-safe
        *
        */
        /*if (!dictFile.exists()) {
            logger.info(String.format("Fasta dictionary file %s does not exist. Trying to create it now.",
                    dictFile.getAbsolutePath()));
            String args[] = {String.format("r=%s", fastaFile.getAbsolutePath()),
                    String.format("o=%s", dictFile.getAbsolutePath())};
            new CreateSequenceDictionary().instanceMain(args);
            logger.info(String.format("Dictionary file created successfully!"));
        }*/

        /**
         * Create fai file if it doesn't exist and throw exception if unsuccessful
         * Note that this implies walkers cannot be run if a fai file is not provided and GATK cannot write to disk
         *
         * Todo: currently working on making this thread safe.
         * Rather than creating a new fai file, structure is created in memory. This whole block will be fixed
         * in a couple days when we figure out the thread stuff.
         *
         */
        /*if (!indexFile.exists()) {
            logger.info(String.format("Fasta index file %s does not exist. Trying to create it now.",
                    indexFile.getAbsolutePath()));
            FastaSequenceIndexBuilder faiBuilder = new FastaSequenceIndexBuilder(fastaFile, this);

            try {
                faiBuilder.saveAsFaiFile();
                logger.info(String.format("Index file created successfully!"));
            }
            catch (Exception e) {
                throw new StingException("Index file could not be saved", e);
            }
        }

        // now create IndexedFastaSequenceFile
        try {
            index = indexFile.exists() ? new IndexedFastaSequenceFile(fastaFile) :
                    new IndexedFastaSequenceFile(fastaFile, sequenceIndex);
        }
        catch (IOException e) {
            throw new StingException("An error occurred while loading the fasta file and its .fasta.fai and " +
                    ".dict counterparts.", e);
        }
        catch (Exception e) {
            throw new StingException("An error occurred while processing the fasta file and its .fasta.fai and " +
                    ".dict counterparts. If the error could have been caused by the .fasta.fai or .dict files, " +
                    "you can re-create them by removing them from the folder that the fasta file is in and " +
                    "running GATK again.", e);
        }*/
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
