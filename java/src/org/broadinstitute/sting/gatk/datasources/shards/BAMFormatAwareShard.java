package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.*;
import net.sf.picard.filter.SamRecordFilter;

import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

/**
 * A common interface for shards that natively understand the BAM format.
 *
 * @author mhanna
 * @version 0.1
 */
public interface BAMFormatAwareShard extends Shard {
    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public Map<SAMReaderID,SAMFileSpan> getFileSpans();

    /**
     * Returns true if this shard is meant to buffer reads, rather
     * than just holding pointers to their locations.
     * @return True if this shard can buffer reads.  False otherwise.
     */
    public boolean buffersReads();

    /**
     * Checks to see whether the buffer is empty.
     * @return True if the buffer is empty.
     */
    public boolean isBufferEmpty();

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    public boolean isBufferFull();

    /**
     * Adds a read to the read buffer.
     * @param read Add a read to the internal shard buffer.
     */
    public void addRead(SAMRecord read);

    /**
     * Assuming this iterator buffers reads, an iterator to the reads
     * stored in the shard.
     * @return An iterator over the reads stored in the shard.
     */
    public StingSAMIterator iterator();

    /**
     * Gets any filter associated with this shard.  Useful for filtering out overlaps, etc.
     * @return A filter if one exists.  Null if not.
     */
    public SamRecordFilter getFilter();
}
