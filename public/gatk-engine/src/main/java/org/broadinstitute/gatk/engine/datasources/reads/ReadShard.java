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
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIteratorAdapter;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

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
     * Default read shard buffer size
     */
    public static final int DEFAULT_MAX_READS = 10000;

    /**
     * What is the maximum number of reads per BAM file which should go into a read shard.
     *
     * TODO: this non-final static variable should either be made final or turned into an
     * TODO: instance variable somewhere -- as both static and mutable it wreaks havoc
     * TODO: with tests that use multiple instances of SAMDataSource (since SAMDataSource
     * TODO: changes this value)
     */
    public static int MAX_READS = DEFAULT_MAX_READS;

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
     *
     * TODO: this mutable static interface is awful and breaks tests -- need to refactor
     *
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
     * Fills this shard's buffer with reads from the iterator passed in
     *
     * @param readIter Iterator from which to draw the reads to fill the shard
     */
    @Override
    public void fill( PeekableIterator<SAMRecord> readIter ) {
        if( ! buffersReads() )
            throw new ReviewedGATKException("Attempting to fill a non-buffering shard.");

        SAMFileHeader.SortOrder sortOrder = getReadProperties().getSortOrder();
        SAMRecord read = null;

        while( ! isBufferFull() && readIter.hasNext() ) {
            final SAMRecord nextRead = readIter.peek();
            if ( read == null || (nextRead.getReferenceIndex().equals(read.getReferenceIndex())) ) {
                // only add reads to the shard if they are on the same contig
                read = readIter.next();
                addRead(read);
            } else {
                break;
            }
        }

        // If the reads are sorted in coordinate order, ensure that all reads
        // having the same alignment start become part of the same shard, to allow
        // downsampling to work better across shard boundaries. Note that because our
        // read stream has already been fed through the positional downsampler, which
        // ensures that at each alignment start position there are no more than dcov
        // reads, we're in no danger of accidentally creating a disproportionately huge
        // shard
        if ( sortOrder == SAMFileHeader.SortOrder.coordinate ) {
            while ( readIter.hasNext() ) {
                SAMRecord additionalRead = readIter.peek();

                // Stop filling the shard as soon as we encounter a read having a different
                // alignment start or contig from the last read added in the earlier loop
                // above, or an unmapped read
                if ( read == null ||
                     additionalRead.getReadUnmappedFlag() ||
                     ! additionalRead.getReferenceIndex().equals(read.getReferenceIndex()) ||
                     additionalRead.getAlignmentStart() != read.getAlignmentStart() ) {
                    break;
                }

                addRead(readIter.next());
            }
        }

        // If the reads are sorted in queryname order, ensure that all reads
        // having the same queryname become part of the same shard.
        if( sortOrder == SAMFileHeader.SortOrder.queryname ) {
            while( readIter.hasNext() ) {
                SAMRecord nextRead = readIter.peek();
                if( read == null || ! read.getReadName().equals(nextRead.getReadName()) )
                    break;
                addRead(readIter.next());
            }
        }
    }

    /**
     * Creates an iterator over reads stored in this shard's read cache.
     * @return
     */
    public GATKSAMIterator iterator() {
        return GATKSAMIteratorAdapter.adapt(reads.iterator());
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
            boolean foundMapped = false;

            for ( final SAMRecord read : reads ) {
                if ( contig != null && ! read.getReferenceName().equals(contig) )
                    throw new ReviewedGATKException("ReadShard contains reads spanning contig boundaries, which is no longer allowed. "
                            + "First contig is " + contig + " next read was " + read.getReferenceName() );
                contig = read.getReferenceName();

                // Even if this shard as a *whole* is not "unmapped", we can still encounter *individual* unmapped mates
                // of mapped reads within this shard's buffer. In fact, if we're very unlucky with shard boundaries,
                // this shard might consist *only* of unmapped mates! We need to refrain from using the alignment
                // starts/stops of these unmapped mates, and detect the case where the shard has been filled *only*
                // with unmapped mates.
                if ( ! read.getReadUnmappedFlag() ) {
                    foundMapped = true;
                    if ( read.getAlignmentStart() < start ) start = read.getAlignmentStart();
                    if ( read.getAlignmentEnd() > stop ) stop = read.getAlignmentEnd();
                }
            }

            assert contig != null;

            if ( ! foundMapped || contig.equals("*") ) // all reads are unmapped
                return GenomeLoc.UNMAPPED;
            else
                return parser.createGenomeLoc(contig, start, stop);
        }
    }
}
