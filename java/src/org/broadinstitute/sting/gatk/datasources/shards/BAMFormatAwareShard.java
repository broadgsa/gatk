package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.*;
import net.sf.picard.filter.SamRecordFilter;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

/**
 * A common interface for shards that natively understand the BAM format.
 *
 * @author mhanna
 * @version 0.1
 */
public abstract class BAMFormatAwareShard extends Shard {
    /**
     * Whether the current location is unmapped.
     */
    private final boolean isUnmapped;

    /**
     * Reads data, if applicable.
     */
    private final SAMDataSource readsDataSource;

    /**
     * The data backing the next chunks to deliver to the traversal engine.
     */
    private final Map<SAMReaderID,SAMFileSpan> fileSpans;

    public BAMFormatAwareShard(GenomeLocParser parser,
                               ShardType shardType,
                               List<GenomeLoc> locs,
                               SAMDataSource readsDataSource,
                               Map<SAMReaderID,SAMFileSpan> fileSpans,
                               boolean isUnmapped) {
        super(parser, shardType, locs);
        this.readsDataSource = readsDataSource;
        this.fileSpans = fileSpans;
        this.isUnmapped = isUnmapped;
    }

    /**
     * Closes the shard, tallying and incorporating read data.
     */
    @Override
    public void close() {
        readsDataSource.incorporateReadMetrics(readMetrics);
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public Map<SAMReaderID,SAMFileSpan> getFileSpans() {
        return Collections.unmodifiableMap(fileSpans);
    }

    /**
     * Gets key read validation and filtering properties.
     * @return set of read properties associated with this shard.
     */
    @Override
    public ReadProperties getReadProperties() {
        return readsDataSource.getReadsInfo();
    }

    /**
     * Returns true if this shard is meant to buffer reads, rather
     * than just holding pointers to their locations.
     * @return True if this shard can buffer reads.  False otherwise.
     */
    public boolean buffersReads() { return false; }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    public boolean isBufferEmpty() { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Returns true if the read buffer is currently full.
     * @return True if this shard's buffer is full (and the shard can buffer reads).
     */
    public boolean isBufferFull() { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Adds a read to the read buffer.
     * @param read Add a read to the internal shard buffer.
     */
    public void addRead(SAMRecord read) { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Gets the iterator over the elements cached in the shard.
     * @return
     */
    public StingSAMIterator iterator() { throw new UnsupportedOperationException("This shard does not buffer reads."); }


    /**
     * Whether this shard points to an unmapped region.
     * Some shard types conceptually be unmapped (e.g. LocusShards).  In
     * this case, isUnmapped should always return false.
     * @return True if this shard is unmapped.  False otherwise.
     */
    public boolean isUnmapped() {
        return isUnmapped;
    }
}
