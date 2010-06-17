package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;

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
     * A list of the chunks associated with this shard.
     */
    private final Map<SAMReaderID,SAMFileSpan> fileSpans;

    // currently our location
    private final List<GenomeLoc> loci;

    /**
     * Create a new locus shard, divided by index.
     * @param intervals List of intervals to process.
     * @param fileSpans File spans associated with that interval.
     */
    public LocusShard(List<GenomeLoc> intervals, Map<SAMReaderID,SAMFileSpan> fileSpans) {
        this.loci = intervals;    
        this.fileSpans = fileSpans;
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
        return ShardType.LOCUS;
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
