package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

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
     * What is the maximum number of reads per BAM file which should go into a read shard.
     */
    public static int MAX_READS = 10000;

    /**
     * The reads making up this shard.
     */
    private final Collection<SAMRecord> reads = new ArrayList<SAMRecord>(MAX_READS);

    public ReadShard(GenomeLocParser parser, SAMDataSource readsDataSource, Map<SAMReaderID,SAMFileSpan> fileSpans, List<GenomeLoc> loci, boolean isUnmapped) {
        super(parser, ShardType.READ, loci, readsDataSource, fileSpans, isUnmapped);
    }

    /**
     * Sets the maximum number of reads buffered in a read shard.  Implemented as a weirdly static interface
     * until we know what effect tuning this parameter has.
     * @param bufferSize New maximum number
     */
    static void setReadBufferSize(final int bufferSize) {
        MAX_READS = bufferSize;
    }

    /**
     * What read buffer size are we using?
     *
     * @return
     */
    public static int getReadBufferSize() {
        return MAX_READS;
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
        return reads.size() > ReadShard.MAX_READS;
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

    /**
     * Get the full span from the start of the left most read to the end of the right most one
     *
     * Note this may be different than the getLocation() of the shard, as this reflects the
     * targeted span, not the actual span of reads
     *
     * @return the genome loc representing the span of these reads on the genome
     */
    public GenomeLoc getReadsSpan() {
        if ( isUnmapped() || super.getGenomeLocs() == null || reads.isEmpty() )
            return super.getLocation();
        else {
            int start = Integer.MAX_VALUE;
            int stop = Integer.MIN_VALUE;
            String contig = null;

            for ( final SAMRecord read : reads ) {
                if ( contig != null && ! read.getReferenceName().equals(contig) )
                    throw new ReviewedStingException("ReadShard contains reads spanning contig boundaries, which is no longer allowed. "
                            + "First contig is " + contig + " next read was " + read.getReferenceName() );
                contig = read.getReferenceName();
                if ( read.getAlignmentStart() < start ) start = read.getAlignmentStart();
                if ( read.getAlignmentEnd() > stop ) stop = read.getAlignmentEnd();
            }

            return parser.createGenomeLoc(contig, start, stop);
        }
    }
}
