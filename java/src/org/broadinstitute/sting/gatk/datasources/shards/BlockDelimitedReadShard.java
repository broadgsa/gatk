package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.*;
import net.sf.picard.filter.SamRecordFilter;

import java.util.*;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

/**
 * Expresses a shard of read data in block format.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDelimitedReadShard extends ReadShard implements BAMFormatAwareShard {
    /**
     * Information about the origins of reads.
     */
    private final Reads sourceInfo;

    /**
     * The data backing the next chunks to deliver to the traversal engine.
     */
    private final Map<SAMReaderID,SAMFileSpan> fileSpans;

    /**
     * The reads making up this shard.
     */
    private final Collection<SAMRecord> reads = new ArrayList<SAMRecord>(BlockDelimitedReadShardStrategy.MAX_READS);

    /**
     * The filter to be applied to all reads meeting this criteria.
     */
    private final SamRecordFilter filter;

    /**
     * An BlockDelimitedLocusShard can be used either for READ or READ shard types.
     * Track which type is being used.
     */
    private final Shard.ShardType shardType;        

    public BlockDelimitedReadShard(Reads sourceInfo, Map<SAMReaderID,SAMFileSpan> fileSpans, SamRecordFilter filter, Shard.ShardType shardType) {
        this.sourceInfo = sourceInfo;
        this.fileSpans = fileSpans;
        this.filter = filter;
        this.shardType = shardType;
    }

    /**
     * Returns true if this shard is meant to buffer reads, rather
     * than just holding pointers to their locations.
     * @return True if this shard can buffer reads.  False otherwise.
     */
    @Override
    public boolean buffersReads() {
        return true;
    }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    @Override
    public boolean isBufferEmpty() {
        return reads.size() == 0;
    }    

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    @Override
    public boolean isBufferFull() {
        return reads.size() > BlockDelimitedReadShardStrategy.MAX_READS; 
    }

    /**
     * Adds a read to the read buffer.
     * @param read Add a read to the internal shard buffer.
     */
    @Override
    public void addRead(SAMRecord read) {
        // DO NOT validate that the buffer is full.  Paired read sharding will occasionally have to stuff another
        // read or two into the buffer.
        reads.add(read);
    }

    /**
     * Creates an iterator over reads stored in this shard's read cache.
     * @return
     */
    @Override
    public StingSAMIterator iterator() {
        return StingSAMIteratorAdapter.adapt(sourceInfo,reads.iterator());
    }

    @Override
    public SamRecordFilter getFilter() {
        return filter;
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    @Override
    public Map<SAMReaderID,SAMFileSpan> getFileSpans() {
        return Collections.unmodifiableMap(fileSpans);
    }

    /**
     * returns the type of shard.
     */
    @Override
    public ShardType getShardType() {
        return shardType;
    }    

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for(Map.Entry<SAMReaderID,SAMFileSpan> entry: fileSpans.entrySet()) {
            sb.append(entry.getKey());
            sb.append(": ");
            sb.append(entry.getValue());
            sb.append(' ');
        }
        return sb.toString();
    }
}
