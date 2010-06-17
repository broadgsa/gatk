package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

import java.util.*;

import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import net.sf.picard.filter.SamRecordFilter;

/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 5:03:13 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */

/**
 * Expresses a shard of read data in block format.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadShard implements BAMFormatAwareShard {
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
    private final Collection<SAMRecord> reads = new ArrayList<SAMRecord>(ReadShardStrategy.MAX_READS);

    /**
     * The filter to be applied to all reads meeting this criteria.
     */
    private final SamRecordFilter filter;

    public ReadShard(Reads sourceInfo, Map<SAMReaderID,SAMFileSpan> fileSpans, SamRecordFilter filter) {
        this.sourceInfo = sourceInfo;
        this.fileSpans = fileSpans;
        this.filter = filter;
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    @Override
    public Map<SAMReaderID,SAMFileSpan> getFileSpans() {
        return Collections.unmodifiableMap(fileSpans);
    }

    /** @return the genome location represented by this shard */
    @Override
    public List<GenomeLoc> getGenomeLocs() {
        throw new UnsupportedOperationException("ReadShard isn't genome loc aware");
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
        return reads.size() > ReadShardStrategy.MAX_READS;
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
     * what kind of shard do we return
     *
     * @return ShardType, indicating the type
     */
    @Override
    public ShardType getShardType() {
        return ShardType.READ;
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
