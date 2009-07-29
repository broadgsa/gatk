/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.executive;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 12:37:23 PM
 * To change this template use File | Settings | File Templates.
 */

/** Shards and schedules data in manageable chunks. */
public abstract class MicroScheduler {
    protected static Logger logger = Logger.getLogger(MicroScheduler.class);

    protected final TraversalEngine traversalEngine;
    protected final IndexedFastaSequenceFile reference;

    private final SAMDataSource reads;
    private final Collection<ReferenceOrderedDataSource> rods;

    /**
     * MicroScheduler factory function.  Create a microscheduler appropriate for reducing the
     * selected walker.
     *
     * @param walker        Which walker to use.
     * @param reads         the informations associated with the reads
     * @param reference     the reference file
     * @param rods          the rods to include in the traversal
     * @param nThreadsToUse Number of threads to utilize.
     *
     * @return The best-fit microscheduler.
     */
    public static MicroScheduler create(Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods, int nThreadsToUse) {
        if (walker instanceof TreeReducible && nThreadsToUse > 1) {
            logger.info("Creating hierarchical microscheduler");
            return new HierarchicalMicroScheduler(walker, reads, reference, rods, nThreadsToUse);
        } else {
            logger.info("Creating linear microscheduler");
            return new LinearMicroScheduler(walker, reads, reference, rods);
        }
    }

    /**
     * Create a microscheduler given the reads and reference.
     *
     * @param walker  the walker to execute with
     * @param reads   The reads.
     * @param reference The reference.
     * @param rods    the rods to include in the traversal
     */
    protected MicroScheduler(Walker walker, SAMDataSource reads, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        if (walker instanceof ReadWalker) {
            traversalEngine = new TraverseReads();
        } else if (walker instanceof LocusWalker) {
            traversalEngine = new TraverseLoci();
        } else if (walker instanceof LocusWindowWalker) {
            traversalEngine = new TraverseLocusWindows();
        } else if (walker instanceof DuplicateWalker) {
            traversalEngine = new TraverseDuplicates();
        } else {
            throw new UnsupportedOperationException("Unable to determine traversal type, the walker is an unknown type.");
        }
        this.reads = reads;
        this.reference = reference;
        this.rods = rods;

        validate(this.reads,this.reference);

        // Side effect: initialize the traversal engine with reads data.
        // TODO: Give users a dedicated way of getting the header so that the MicroScheduler
        //       doesn't have to bend over backward providing legacy getters and setters.
        traversalEngine.setSAMHeader(reads.getHeader());        
        traversalEngine.initialize();                
    }

    /**
     * A temporary getter for the traversal engine.  In the future, clients
     * of the microscheduler shouldn't need to know anything about the traversal engine.
     *
     * @return The traversal engine.
     */
    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }    

    /**
     * Walks a walker over the given list of intervals.
     *
     * @param walker        Computation to perform over dataset.
     * @param shardStrategy A strategy for sharding the data.
     *
     * @return the return type of the walker
     */
    public abstract Object execute(Walker walker, ShardStrategy shardStrategy);


    /**
     * Gets an window into all the data that can be viewed as a single shard.
     *
     * @param shard The section of data to view.
     *
     * @return An accessor for all the data in this shard.
     */
    protected ShardDataProvider getShardDataProvider(Shard shard) {
        return new ShardDataProvider(shard, reads, reference, rods);
    }

    /**
     * Print summary information for the analysis.
     * @param sum The final reduce output.
     */
    protected void printOnTraversalDone(Object sum) {
        // HACK: The microscheduler should be too dumb to know anything about the data
        //       it's actually processing; it should just funnel anything it receives
        //       to the traversal engine.
        // TODO: Implement code to allow the datasources to print summary info of the
        //       data they've seen.
        if( reads != null && reads.getViolationHistogram().getViolationCount() > 0 )
            logger.warn(String.format("%n%s",reads.getViolationHistogram()));

        traversalEngine.printOnTraversalDone(sum);
    }

    /**
     * Returns data source maintained by this scheduler
     * @return
     */
    public SAMDataSource getSAMDataSource() { return reads; }

    /**
     * Returns the reference maintained by this scheduler.
     * @return The reference maintained by this scheduler.
     */
    public IndexedFastaSequenceFile getReference() { return reference; }

    /**
     * Now that all files are open, validate the sequence dictionaries of the reads vs. the reference.
     * TODO: Doing this in the MicroScheduler is a bit late, but this is where data sources are initialized.
     * TODO: Move the initialization of data sources back to the GenomeAnalysisEngine.
     * @param reads Reads data source.
     * @param reference Reference data source.
     */
    private void validate( SAMDataSource reads, ReferenceSequenceFile reference ) {
        if( reads == null || reference == null )
            return;
        
        // Compile a set of sequence names that exist in the BAM files.
        SAMSequenceDictionary readsDictionary = reads.getHeader().getSequenceDictionary();

        Set<String> readsSequenceNames = new TreeSet<String>();
        for( SAMSequenceRecord dictionaryEntry: readsDictionary.getSequences() )
            readsSequenceNames.add(dictionaryEntry.getSequenceName());

        // Compile a set of sequence names that exist in the reference file.
        SAMSequenceDictionary referenceDictionary = reference.getSequenceDictionary();

        Set<String> referenceSequenceNames = new TreeSet<String>();
        for( SAMSequenceRecord dictionaryEntry: referenceDictionary.getSequences() )
            referenceSequenceNames.add(dictionaryEntry.getSequenceName());

        if( readsSequenceNames.size() == 0 ) {
            logger.info("Reads file is unmapped.  Skipping validation against reference.");
            return;
        }

        // If there's no overlap between reads and reference, data will be bogus.  Throw an exception.
        Set<String> intersectingSequenceNames = new HashSet<String>(readsSequenceNames);
        intersectingSequenceNames.retainAll(referenceSequenceNames);
        if( intersectingSequenceNames.size() == 0 ) {
            StringBuilder error = new StringBuilder();
            error.append("No overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.  Perhaps you're using the wrong reference?\n");
            error.append(System.getProperty("line.separator"));
            error.append(String.format("Reads contigs:     %s%n", prettyPrintSequenceRecords(readsDictionary)));
            error.append(String.format("Reference contigs: %s%n", prettyPrintSequenceRecords(referenceDictionary)));
            logger.error(error.toString());
            Utils.scareUser("No overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.");
        }

        // If the two datasets are not equal and neither is a strict subset of the other, warn the user.
        if( !readsSequenceNames.equals(referenceSequenceNames) &&
            !readsSequenceNames.containsAll(referenceSequenceNames) &&
            !referenceSequenceNames.containsAll(readsSequenceNames)) {
            StringBuilder warning = new StringBuilder();
            warning.append("Limited overlap exists between sequence dictionary of the reads and the sequence dictionary of the reference.  Perhaps you're using the wrong reference?\n");
            warning.append(System.getProperty("line.separator"));
            warning.append(String.format("Reads contigs:     %s%n", prettyPrintSequenceRecords(readsDictionary)));
            warning.append(String.format("Reference contigs: %s%n", prettyPrintSequenceRecords(referenceDictionary)));
            logger.warn(warning.toString());
        }
    }

    private String prettyPrintSequenceRecords( SAMSequenceDictionary sequenceDictionary ) {
        String[] sequenceRecordNames = new String[ sequenceDictionary.size() ];
        int sequenceRecordIndex = 0;
        for( SAMSequenceRecord sequenceRecord: sequenceDictionary.getSequences() )
            sequenceRecordNames[sequenceRecordIndex++] = sequenceRecord.getSequenceName();
        return Arrays.deepToString(sequenceRecordNames);
    }
}
