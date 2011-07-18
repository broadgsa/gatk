package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;

import java.util.Collections;
import java.util.List;
import java.util.Map;
/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 5:00:27 PM
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
 * @author aaron
 * @version 1.0
 * @date Apr 10, 2009
 * <p/>
 * Interface Shard
 * <p/>
 * The base abstract class for shards.
 */
public abstract class Shard implements HasGenomeLocation {
    public enum ShardType {
        READ, LOCUS
    }

    protected final GenomeLocParser parser; // incredibly annoying!

    /**
     * What type of shard is this?  Read or locus?
     */
    protected final ShardType shardType;

    /**
     * Locations.
     */
    protected final List<GenomeLoc> locs;

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



    /**
     * Statistics about which reads in this shards were used and which were filtered away.
     */
    protected final ReadMetrics readMetrics = new ReadMetrics();

    /**
     * Whether this shard points to an unmapped region.
     * Some shard types conceptually be unmapped (e.g. LocusShards).  In
     * this case, isUnmapped should always return false.
     * @return True if this shard is unmapped.  False otherwise.
     */
    public boolean isUnmapped() {
        return isUnmapped;
    }    

    public Shard(GenomeLocParser parser,
                 ShardType shardType,
                 List<GenomeLoc> locs,
                 SAMDataSource readsDataSource,
                 Map<SAMReaderID,SAMFileSpan> fileSpans,
                 boolean isUnmapped) {
        this.locs = locs;
        this.parser = parser;
        this.shardType = shardType;
        this.readsDataSource = readsDataSource;
        this.fileSpans = fileSpans;
        this.isUnmapped = isUnmapped;        
    }

    /**
     * If isUnmapped is true, than getGenomeLocs by
     * definition will return a singleton list with a GenomeLoc.UNMAPPED
     *
     * Can return null, indicating that the entire genome is covered.
     *
     * @return the genome location represented by this shard
     */
    public List<GenomeLoc> getGenomeLocs() {
        return locs;
    }

    /**
     * Get the list of chunks delimiting this shard.
     * @return a list of chunks that contain data for this shard.
     */
    public Map<SAMReaderID,SAMFileSpan> getFileSpans() {
        return Collections.unmodifiableMap(fileSpans);
    }    

    /**
     * Returns the span of the genomeLocs comprising this shard
     * @param
     * @return
     */
    public GenomeLoc getLocation() {
        if ( getGenomeLocs() == null )
            return GenomeLoc.WHOLE_GENOME;

        int start = Integer.MAX_VALUE;
        int stop = Integer.MIN_VALUE;
        String contig = null;

        for ( GenomeLoc loc : getGenomeLocs() ) {
            if ( GenomeLoc.isUnmapped(loc) )
                // special case the unmapped region marker, just abort out
                return loc;
            contig = loc.getContig();
            if ( loc.getStart() < start ) start = loc.getStart();
            if ( loc.getStop() > stop ) stop = loc.getStop();
        }

        return parser.createGenomeLoc(contig, start, stop);
    }


    /**
     * what kind of shard do we return
     * @return ShardType, indicating the type
     */
    public ShardType getShardType() {
        return shardType;
    }

    /**
     * Does any releasing / aggregation required when the shard is through being processed.
     */
    public void close() {
        readsDataSource.incorporateReadMetrics(readMetrics);
    }

    /**
     * Gets key read validation and filtering properties.
     * @return set of read properties associated with this shard.
     */
    public ReadProperties getReadProperties() {
        return readsDataSource.getReadsInfo();
    }

    /**
     * Gets the runtime metrics associated with this shard.
     * Retrieves a storage space of metrics about number of reads included, filtered, etc.
     * @return Storage space for metrics.
     */
    public ReadMetrics getReadMetrics() {
        return readMetrics;
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
}
