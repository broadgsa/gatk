package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.BAMFileSpan;
import net.sf.picard.filter.SamRecordFilter;

import java.util.List;
import java.util.Map;


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

/**
 * A shard that's delimited based on the index rather than
 */
public class IndexDelimitedLocusShard extends LocusShard implements BAMFormatAwareShard {
    /**
     * A list of the chunks associated with this shard.
     */
    private final Map<SAMReaderID, BAMFileSpan> fileSpans;

    /**
     * An IndexDelimitedLocusShard can be used either for LOCUS or LOCUS_INTERVAL shard types.
     * Track which type is being used.
     */
    private final ShardType shardType;

    /**
     * Create a new locus shard, divided by index.
     * @param intervals List of intervals to process.
     * @param fileSpans File spans associated with that interval.
     * @param shardType Type of the shard; must be either LOCUS or LOCUS_INTERVAL.
     */
    IndexDelimitedLocusShard(List<GenomeLoc> intervals, Map<SAMReaderID,BAMFileSpan> fileSpans, ShardType shardType) {
        super(intervals);
        this.fileSpans = fileSpans;
        if(shardType != ShardType.LOCUS && shardType != ShardType.LOCUS_INTERVAL)
            throw new StingException("Attempted to create an IndexDelimitedLocusShard with invalid shard type: " + shardType);
        this.shardType = shardType;
    }

    /**
     * Gets the file spans associated with this locus shard.
     * @return A list of the file spans to use when retrieving locus data.
     */
    @Override
    public Map<SAMReaderID,BAMFileSpan> getFileSpans() {
        return fileSpans;
    }

    /**
     * Returns true if this shard is meant to buffer reads, rather
     * than just holding pointers to their locations.
     * @return True if this shard can buffer reads.  False otherwise.
     */
    @Override
    public boolean buffersReads() { return false; }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    @Override
    public boolean isBufferEmpty() { throw new UnsupportedOperationException("This shard does not buffer reads."); }    

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    @Override
    public boolean isBufferFull() { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Adds a read to the read buffer.
     * @param read Add a read to the internal shard buffer.
     */
    @Override
    public void addRead(SAMRecord read) { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Gets the iterator over the elements cached in the shard.
     * @return
     */
    @Override
    public StingSAMIterator iterator() { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Gets a filter testing for overlap of this read with the given shard.
     * @return A filter capable of filtering out reads outside a given shard.
     */
    @Override
    public SamRecordFilter getFilter() {
        return new ReadOverlapFilter(loci);    
    }

    /**
     * returns the type of shard.
     */
    @Override
    public ShardType getShardType() {
        return shardType;
    }
}