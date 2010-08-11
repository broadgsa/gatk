/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.ReadProperties;

import java.util.*;
import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileSpan;
import net.sf.picard.reference.IndexedFastaSequenceFile;

/**
 * A sharding strategy for loci based on reading of the index.
 */
public class LocusShardStrategy implements ShardStrategy {
    /**
     * The data source to use when performing this sharding.
     */
    private final SAMDataSource reads;

    /**
     * An iterator through the available file pointers.
     */
    private final Iterator<FilePointer> filePointerIterator;

    /**
     * construct the shard strategy from a seq dictionary, a shard size, and and genomeLocs
     * @param reads Data source from which to load index data.
     * @param locations List of locations for which to load data.
     */
    LocusShardStrategy(SAMDataSource reads, IndexedFastaSequenceFile reference, GenomeLocSortedSet locations) {
        this.reads = reads;
        if(!reads.isEmpty()) {
            List<GenomeLoc> intervals;
            if(locations == null) {
                // If no locations were passed in, shard the entire BAM file.
                SAMFileHeader header = reads.getHeader();
                intervals = new ArrayList<GenomeLoc>();

                for(SAMSequenceRecord readsSequenceRecord: header.getSequenceDictionary().getSequences()) {
                    // Check this sequence against the reference sequence dictionary.
                    // TODO: Do a better job of merging reads + reference.
                    SAMSequenceRecord refSequenceRecord = reference.getSequenceDictionary().getSequence(readsSequenceRecord.getSequenceName());
                    if(refSequenceRecord != null) {
                        final int length = Math.min(readsSequenceRecord.getSequenceLength(),refSequenceRecord.getSequenceLength());
                        intervals.add(GenomeLocParser.createGenomeLoc(readsSequenceRecord.getSequenceName(),1,length));
                    }
                }
            }
            else
                intervals = locations.toList();

            this.filePointerIterator = IntervalSharder.shardIntervals(this.reads,intervals);
        }
        else {
            final int maxShardSize = 100000;
            List<FilePointer> filePointers = new ArrayList<FilePointer>();
            if(locations == null) {
                for(SAMSequenceRecord refSequenceRecord: reference.getSequenceDictionary().getSequences()) {
                    for(int shardStart = 1; shardStart <= refSequenceRecord.getSequenceLength(); shardStart += maxShardSize) {
                        final int shardStop = Math.min(shardStart+maxShardSize-1, refSequenceRecord.getSequenceLength());
                        filePointers.add(new FilePointer(GenomeLocParser.createGenomeLoc(refSequenceRecord.getSequenceName(),shardStart,shardStop)));
                    }
                }
            }
            else {
                for(GenomeLoc interval: locations)
                    filePointers.add(new FilePointer(interval));
            }
            filePointerIterator = filePointers.iterator();
        }

    }

    /**
     * returns true if there are additional shards
     *
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return filePointerIterator.hasNext();
    }

    /**
     * gets the next Shard
     *
     * @return the next shard
     */
    public LocusShard next() {
        FilePointer nextFilePointer = filePointerIterator.next();
        Map<SAMReaderID,SAMFileSpan> fileSpansBounding = nextFilePointer.fileSpans != null ? nextFilePointer.fileSpans : null;
        return new LocusShard(reads,nextFilePointer.locations,fileSpansBounding);
    }

    /** we don't support the remove command */
    public void remove() {
        throw new UnsupportedOperationException("ShardStrategies don't support remove()");
    }

    /**
     * makes the IntervalShard iterable, i.e. usable in a for loop.
     *
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }
}