package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.ReadProperties;

import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import net.sf.picard.filter.SamRecordFilter;

/**
 * Handles locus shards of BAM information.
 * @author aaron
 * @version 1.0
 * @date Apr 7, 2009
 */
public class LocusShard implements BAMFormatAwareShard {
    /**
     * Source for read data.
     */
    private SAMDataSource dataSource;

    /**
     * A list of the chunks associated with this shard.
     */
    private final Map<SAMReaderID,SAMFileSpan> fileSpans;

    // currently our location
    private final List<GenomeLoc> loci;

    /**
     * Statistics about which reads in this shards were used and which were filtered away.
     */
    private final ReadMetrics readMetrics = new ReadMetrics();

    /**
     * Create a new locus shard, divided by index.
     * @param intervals List of intervals to process.
     * @param fileSpans File spans associated with that interval.
     */
    public LocusShard(SAMDataSource dataSource, List<GenomeLoc> intervals, Map<SAMReaderID,SAMFileSpan> fileSpans) {
        this.dataSource = dataSource;
        this.loci = intervals;    
        this.fileSpans = fileSpans;
    }

    /**
     * Closes the shard, tallying and incorporating read data.
     */
    @Override
    public void close() {
        dataSource.incorporateReadMetrics(readMetrics);
    }

    /**
     * Gets the file spans associated with this locus shard.
     * @return A list of the file spans to use when retrieving locus data.
     */
    @Override
    public Map<SAMReaderID,SAMFileSpan> getFileSpans() {
        return fileSpans;
    }

    /** @return the genome location represented by this shard */
    public List<GenomeLoc> getGenomeLocs() {
        return loci;
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
     * returns the type of shard.
     */
    @Override
    public ShardType getShardType() {
        return ShardType.LOCUS;
    }

    /**
     * Locus shards don't make sense as unmapped regions.  Always return false.
     * @return False always.
     */
    @Override
    public boolean isUnmapped() {
        return false;
    }

    /**
     * Gets key read validation and filtering properties.
     * @return set of read properties associated with this shard.
     */
    @Override
    public ReadProperties getReadProperties() {
        return dataSource.getReadsInfo();
    }

    /**
     * Retrieves a storage space of metrics about number of reads included, filtered, etc.
     * @return Storage space for metrics.
     */    
    @Override
    public ReadMetrics getReadMetrics() {
        return readMetrics;
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return Utils.join(";",loci);
    }
}
