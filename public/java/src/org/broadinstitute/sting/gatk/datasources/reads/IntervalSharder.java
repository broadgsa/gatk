/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.Iterator;

/**
 * Handles the process of aggregating BAM intervals into individual shards.
 * TODO: The task performed by IntervalSharder is now better performed by LocusShardBalancer.  Merge BAMScheduler and IntervalSharder.
 */
public class IntervalSharder implements Iterator<FilePointer> {
    /**
     * The iterator actually laying out the data for BAM scheduling.
     */
    private final PeekableIterator<FilePointer> wrappedIterator;

    /**
     * The parser, for interval manipulation.
     */
    private final GenomeLocParser parser;

    public static IntervalSharder shardOverAllReads(final SAMDataSource dataSource, final GenomeLocParser parser) {
        return new IntervalSharder(BAMScheduler.createOverAllReads(dataSource,parser),parser);
    }

    public static IntervalSharder shardOverMappedReads(final SAMDataSource dataSource, final SAMSequenceDictionary sequenceDictionary, final GenomeLocParser parser) {
        return new IntervalSharder(BAMScheduler.createOverMappedReads(dataSource,sequenceDictionary,parser),parser);
    }

    public static IntervalSharder shardOverIntervals(final SAMDataSource dataSource, final GenomeLocSortedSet loci) {
        return new IntervalSharder(BAMScheduler.createOverIntervals(dataSource,loci),loci.getGenomeLocParser());
    }

    private IntervalSharder(final BAMScheduler scheduler, final GenomeLocParser parser) {
        wrappedIterator = new PeekableIterator<FilePointer>(scheduler);
        this.parser = parser;
    }

    public boolean hasNext() {
        return wrappedIterator.hasNext();
    }

    /**
     * Accumulate shards where there's no additional cost to processing the next shard in the sequence.
     * @return The next file pointer to process.
     */
    public FilePointer next() {
        FilePointer current = wrappedIterator.next();

        while ( wrappedIterator.hasNext() &&
                current.isRegionUnmapped == wrappedIterator.peek().isRegionUnmapped &&
                (current.getContigIndex() == wrappedIterator.peek().getContigIndex() || current.isRegionUnmapped) &&
                current.minus(wrappedIterator.peek()) == 0 ) {

            current = current.combine(parser,wrappedIterator.next());
        }

        return current;
    }

    public void remove() { throw new UnsupportedOperationException("Unable to remove from an interval sharder."); }
}
