/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMFileSpan;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.engine.ReadMetrics;
import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

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
     * Lazy-calculated span of all of the genome locs in this shard
     */
    private GenomeLoc spanningLocation = null;

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
     * @return a GenomeLoc that starts as the first position in getGenomeLocs() and stops at the stop of the last
     *    position in getGenomeLocs()
     */
    public GenomeLoc getLocation() {
        if ( spanningLocation == null ) {
            if ( getGenomeLocs() == null )
                spanningLocation = GenomeLoc.WHOLE_GENOME;
            else if ( getGenomeLocs().size() == 0 ) {
                spanningLocation = getGenomeLocs().get(0);
            } else {
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

                spanningLocation = parser.createGenomeLoc(contig, start, stop);
            }
        }

        return spanningLocation;
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
     * Fills the shard with reads. Can only do this with shards that buffer reads
     * @param readIter Iterator from which to draw the reads to fill the shard
     */
    public void fill( PeekableIterator<SAMRecord> readIter ) { throw new UnsupportedOperationException("This shard does not buffer reads."); }

    /**
     * Gets the iterator over the elements cached in the shard.
     * @return
     */
    public GATKSAMIterator iterator() { throw new UnsupportedOperationException("This shard does not buffer reads."); }    
}
