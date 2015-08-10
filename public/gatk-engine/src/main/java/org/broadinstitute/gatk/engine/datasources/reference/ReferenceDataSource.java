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

package org.broadinstitute.gatk.engine.datasources.reference;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.engine.datasources.reads.LocusShard;
import org.broadinstitute.gatk.engine.datasources.reads.SAMDataSource;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

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
        reference = CachingIndexedFastaSequenceFile.checkAndCreate(fastaFile);
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
