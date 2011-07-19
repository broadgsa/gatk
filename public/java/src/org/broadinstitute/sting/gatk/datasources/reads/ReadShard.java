package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

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
public class ReadShard extends Shard {
    /**
     * The reads making up this shard.
     */
    private final Collection<SAMRecord> reads = new ArrayList<SAMRecord>(ReadShardStrategy.MAX_READS);

    public ReadShard(GenomeLocParser parser, SAMDataSource readsDataSource, Map<SAMReaderID,SAMFileSpan> fileSpans, List<GenomeLoc> loci, boolean isUnmapped) {
        super(parser, ShardType.READ, loci, readsDataSource, fileSpans, isUnmapped);
    }

    /**
     * Returns true if this shard is meant to buffer reads, rather
     * than just holding pointers to their locations.
     * @return True if this shard can buffer reads.  False otherwise.
     */
    public boolean buffersReads() {
        return true;
    }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    public boolean isBufferEmpty() {
        return reads.size() == 0;
    }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    public boolean isBufferFull() {
        return reads.size() > ReadShardStrategy.MAX_READS;
    }

    /**
     * Adds a read to the read buffer.
     * @param read Add a read to the internal shard buffer.
     */
    public void addRead(SAMRecord read) {
        // DO NOT validate that the buffer is full.  Paired read sharding will occasionally have to stuff another
        // read or two into the buffer.
        reads.add(read);
    }

    /**
     * Creates an iterator over reads stored in this shard's read cache.
     * @return
     */
    public StingSAMIterator iterator() {
        return StingSAMIteratorAdapter.adapt(reads.iterator());
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for(Map.Entry<SAMReaderID,SAMFileSpan> entry: getFileSpans().entrySet()) {
            sb.append(entry.getKey());
            sb.append(": ");
            sb.append(entry.getValue());
            sb.append(' ');
        }
        return sb.toString();
    }
}
