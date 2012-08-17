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
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.gatk.datasources.reads.LocusShard;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.file.FSLockWithShared;
import org.broadinstitute.sting.utils.file.FileSystemInabilityToLockException;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Loads reference data from fasta file
 * Looks for fai and dict files, and tries to create them if they don't exist
 */
public class ReferenceDataSource {
    private IndexedFastaSequenceFile reference;

    /** our log, which we want to capture anything from this class */
    protected static final org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(ReferenceDataSource.class);

    /**
     * Create reference data source from fasta file
     * @param fastaFile Fasta file to be used as reference
     */
    public ReferenceDataSource(File fastaFile) {
        // does the fasta file exist? check that first...
        if (!fastaFile.exists())
            throw new UserException("The fasta file you specified (" + fastaFile.getAbsolutePath() + ") does not exist.");

        final boolean isGzipped = fastaFile.getAbsolutePath().endsWith(".gz");

        final File indexFile = new File(fastaFile.getAbsolutePath() + ".fai");

        // determine the name for the dict file
        final String fastaExt = (fastaFile.getAbsolutePath().endsWith("fa") ? ".fa" : ".fasta" ) + (isGzipped ? ".gz" : "");
        final File dictFile = new File(fastaFile.getAbsolutePath().replace(fastaExt, ".dict"));

        /*
        * if index file does not exist, create it manually
        */
        if (!indexFile.exists()) {
            if ( isGzipped ) throw new UserException.CouldNotCreateReferenceFAIorDictForGzippedRef(fastaFile);

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
            if ( isGzipped ) throw new UserException.CouldNotCreateReferenceFAIorDictForGzippedRef(fastaFile);

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

            reference = new CachingIndexedFastaSequenceFile(fastaFile);

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
        return this.reference;
    }

    /**
     * Creates an iterator for processing the entire reference.
     * @param readsDataSource the reads datasource to embed in the locus shard.
     * @param parser used to generate/regenerate intervals.  TODO: decouple the creation of the shards themselves from the creation of the driving iterator so that datasources need not be passed to datasources.
     * @param maxShardSize The maximum shard size which can be used to create this list.
     * @return Creates a schedule for performing a traversal over the entire reference.
     */
    public Iterable<Shard> createShardsOverEntireReference(final SAMDataSource readsDataSource, final GenomeLocParser parser, final int maxShardSize) {
        List<Shard> shards = new ArrayList<Shard>();
        for(SAMSequenceRecord refSequenceRecord: reference.getSequenceDictionary().getSequences()) {
            for(int shardStart = 1; shardStart <= refSequenceRecord.getSequenceLength(); shardStart += maxShardSize) {
                final int shardStop = Math.min(shardStart+maxShardSize-1, refSequenceRecord.getSequenceLength());
                shards.add(new LocusShard(parser,
                        readsDataSource,
                        Collections.singletonList(parser.createGenomeLoc(refSequenceRecord.getSequenceName(),shardStart,shardStop)),
                        null));
            }
        }
        return shards;
    }


    public Iterable<Shard> createShardsOverIntervals(final SAMDataSource readsDataSource, final GenomeLocSortedSet intervals, final int maxShardSize) {
        List<Shard> shards = new ArrayList<Shard>();

        for(GenomeLoc interval: intervals) {
            while(interval.size() > maxShardSize) {
                shards.add(new LocusShard(intervals.getGenomeLocParser(),
                        readsDataSource,
                        Collections.singletonList(intervals.getGenomeLocParser().createGenomeLoc(interval.getContig(),interval.getStart(),interval.getStart()+maxShardSize-1)),
                        null));
                interval = intervals.getGenomeLocParser().createGenomeLoc(interval.getContig(),interval.getStart()+maxShardSize,interval.getStop());
            }
            shards.add(new LocusShard(intervals.getGenomeLocParser(),
                    readsDataSource,
                    Collections.singletonList(interval),
                    null));
        }

        return shards;
    }


    /**
     * Creates an iterator for processing the entire reference.
     * @param readsDataSource  the reads datasource to embed in the locus shard.  TODO: decouple the creation of the shards themselves from the creation of the driving iterator so that datasources need not be passed to datasources.
     * @param intervals        the list of intervals to use when processing the reference.
     * @param targetShardSize  the suggested - and maximum - shard size which can be used to create this list; we will merge intervals greedily so that we generate shards up to but not greater than the target size.
     * @return Creates a schedule for performing a traversal over the entire reference.
     */
/*
    public Iterable<Shard> createShardsOverIntervals(final SAMDataSource readsDataSource, final GenomeLocSortedSet intervals, final int targetShardSize) {
        final List<Shard> shards = new ArrayList<Shard>();
        final GenomeLocParser parser = intervals.getGenomeLocParser();
        LinkedList<GenomeLoc> currentIntervals = new LinkedList<GenomeLoc>();

        for(GenomeLoc interval: intervals) {
            // if the next interval is too big, we can safely shard currentInterval and then break down this one
            if (interval.size() > targetShardSize) {
                if (!currentIntervals.isEmpty())
                    shards.add(createShardFromInterval(currentIntervals, readsDataSource, parser));
                while(interval.size() > targetShardSize) {
                    final GenomeLoc partialInterval = parser.createGenomeLoc(interval.getContig(), interval.getStart(), interval.getStart()+targetShardSize-1);
                    shards.add(createShardFromInterval(Collections.singletonList(partialInterval), readsDataSource, parser));
                    interval = parser.createGenomeLoc(interval.getContig(), interval.getStart() + targetShardSize, interval.getStop());
                }
                currentIntervals = new LinkedList<GenomeLoc>();
                currentIntervals.add(interval);
            }
            // otherwise, we need to check whether we can merge this interval with currentInterval (and either shard currentInterval or merge accordingly)
            else {
                if (currentIntervals.isEmpty()) {
                    currentIntervals.add(interval);
                }
                else {
                    if (currentIntervals.getLast().compareContigs(interval) != 0 || interval.getStop() - currentIntervals.getLast().getStart() + 1 > targetShardSize) {
                        shards.add(createShardFromInterval(currentIntervals, readsDataSource, parser));
                        currentIntervals = new LinkedList<GenomeLoc>();
                    }
                    currentIntervals.add(interval);
                }
            }
        }
        if (!currentIntervals.isEmpty())
            shards.add(createShardFromInterval(currentIntervals, readsDataSource, parser));
        return shards;
    }

    private static Shard createShardFromInterval(final List<GenomeLoc> intervals, final SAMDataSource readsDataSource, final GenomeLocParser parser) {
        //logger.debug("Adding shard " + interval);
        return new LocusShard(parser,
                readsDataSource,
                intervals,
                null);
    }
*/
}
