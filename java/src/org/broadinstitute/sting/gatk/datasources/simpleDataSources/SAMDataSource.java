package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.shards.ReadShard;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.SAMReadValidator;
import org.broadinstitute.sting.utils.sam.SAMReadViolationHistogram;

import java.io.File;
import java.util.List;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * User: aaron
 * Date: Mar 26, 2009
 * Time: 2:36:16 PM
 * <p/>
 * Converts shards to SAM iterators over the specified region
 */
public class SAMDataSource implements SimpleDataSource {


    /** Backing support for reads. */
    private final Reads reads;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMDataSource.class);

    // used for the reads case, the last count of reads retrieved
    long readsTaken = 0;

    // our last genome loc position
    GenomeLoc lastReadPos = null;

    // do we take unmapped reads
    private boolean includeUnmappedReads = true;

    // reads based traversal variables
    private boolean intoUnmappedReads = false;
    private int readsSeenAtLastPos = 0;

    // A pool of SAM iterators.
    private SAMIteratorPool iteratorPool = null;

    /**
     * Returns a histogram of reads that were screened out, grouped by the nature of the error.
     * @return Histogram of reads.  Will not be null.
     */
    public SAMReadViolationHistogram getViolationHistogram() {
        return iteratorPool.getViolationHistogram();
    }

    /**
     * constructor, given sam files
     *
     * @param reads the list of sam files
     */
    public SAMDataSource( Reads reads ) throws SimpleDataSourceLoadException {
        this.reads = reads;

        // check the length
        if (reads.getReadsFiles().size() < 1) {
            throw new SimpleDataSourceLoadException("SAMDataSource: you must provide a list of length greater then 0");
        }
        for (File smFile : reads.getReadsFiles()) {
            if (!smFile.canRead()) {
                throw new SimpleDataSourceLoadException("SAMDataSource: Unable to load file: " + smFile.getName());
            }
        }
        iteratorPool = new SAMIteratorPool(reads);
    }

    /**
     * Gets the (potentially merged) SAM file header.
     *
     * @return SAM file header.
     */
    public SAMFileHeader getHeader() {
        return iteratorPool.getHeader();
    }

    /**
     * 
     * @param shard the shard to get data for
     *
     * @return an iterator for that region
     */
    public StingSAMIterator seek( Shard shard ) throws SimpleDataSourceLoadException {
        // setup the iterator pool if it's not setup
        boolean queryOverlapping = ( shard.getShardType() == Shard.ShardType.READ ) ? false : true;
        iteratorPool.setQueryOverlapping(queryOverlapping);

        StingSAMIterator iterator = null;
        if (shard.getShardType() == Shard.ShardType.READ) {
            iterator = seekRead((ReadShard) shard);
            iterator = applyDecoratingIterators(true,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getFilterZeroMappingQualityReads(),
                    reads.getSafetyChecking());
        } else if (shard.getShardType() == Shard.ShardType.LOCUS || shard.getShardType() == Shard.ShardType.INTERVAL) {
            iterator = seekLocus(shard.getGenomeLoc());
            iterator = applyDecoratingIterators(false,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getFilterZeroMappingQualityReads(),
                    reads.getSafetyChecking());
        } else {
            throw new StingException("seek: Unknown shard type");
        }

        return iterator;
    }


    /**
     * <p>
     * seekLocus
     * </p>
     *
     * @param location the genome location to extract data for
     *
     * @return an iterator for that region
     */
    private StingSAMIterator seekLocus( GenomeLoc location ) throws SimpleDataSourceLoadException {
        return iteratorPool.iterator(new MappedStreamSegment(location));
    }


    /**
     * <p>
     * seek
     * </p>
     *
     * @param shard the read shard to extract from
     *
     * @return an iterator for that region
     */
    private StingSAMIterator seekRead( ReadShard shard ) throws SimpleDataSourceLoadException {
        StingSAMIterator iter = null;

        if (!intoUnmappedReads) {
            if (lastReadPos == null) {
                lastReadPos = GenomeLocParser.createGenomeLoc(getHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, Integer.MAX_VALUE);
                iter = iteratorPool.iterator(new MappedStreamSegment(lastReadPos));
                return InitialReadIterator(shard.getSize(), iter);
            } else {
                lastReadPos = GenomeLocParser.setStop(lastReadPos,-1);
                iter = fastMappedReadSeek(shard.getSize(), StingSAMIteratorAdapter.adapt(reads, iteratorPool.iterator(new MappedStreamSegment(lastReadPos))));
            }

            if (intoUnmappedReads && !includeUnmappedReads)
                shard.signalDone();
        }

        if (intoUnmappedReads && includeUnmappedReads) {
            if (iter != null)
                iter.close();
            iter = toUnmappedReads(shard.getSize());
            if (!iter.hasNext())
                shard.signalDone();
        }

        return iter;
    }

    /**
     * If we're in by-read mode, this indicates if we want
     * to see unmapped reads too.  Only seeing mapped reads
     * is much faster, but most BAM files have significant
     * unmapped read counts.
     *
     * @param seeUnMappedReads true to see unmapped reads, false otherwise
     */
    public void viewUnmappedReads( boolean seeUnMappedReads ) {
        includeUnmappedReads = seeUnMappedReads;
    }

    /**
     * For unit testing, add a custom iterator pool.
     *
     * @param iteratorPool Custom mock iterator pool.
     */
    void setResourcePool( SAMIteratorPool iteratorPool ) {
        this.iteratorPool = iteratorPool;
    }    

    /**
     * Retrieve unmapped reads.
     *
     * @param readCount how many reads to retrieve
     *
     * @return the bounded iterator that you can use to get the intervaled reads from
     */
    StingSAMIterator toUnmappedReads( long readCount ) {
        StingSAMIterator iter = iteratorPool.iterator(new UnmappedStreamSegment(readsTaken, readCount));
        readsTaken += readCount;
        return iter;
    }


    /**
     * A seek function for mapped reads.
     *
     * @param readCount how many reads to retrieve
     * @param iter      the iterator to use, seeked to the correct start location
     *
     * @return the bounded iterator that you can use to get the intervaled reads from.  Will be a zero-length
     *         iterator if no reads are available.
     * @throws SimpleDataSourceLoadException
     */
    StingSAMIterator fastMappedReadSeek( long readCount, StingSAMIterator iter ) throws SimpleDataSourceLoadException {
        BoundedReadIterator bound;
        correctForReadPileupSeek(iter);
        if (readsTaken == 0) {
            return InitialReadIterator(readCount, iter);
        }
        int x = 0;
        SAMRecord rec = null;
        int lastPos = 0;

        while (x < readsTaken) {
            if (iter.hasNext()) {
                rec = iter.next();
                if (lastPos == rec.getAlignmentStart()) ++this.readsSeenAtLastPos;
                else this.readsSeenAtLastPos = 1;
                lastPos = rec.getAlignmentStart();
                ++x;
            } else {
                iter.close();

                // jump contigs
                lastReadPos = GenomeLocParser.toNextContig(lastReadPos);
                if (lastReadPos == null) {
                    // check to see if we're using unmapped reads, if not return, we're done
                    readsTaken = 0;
                    intoUnmappedReads = true;

                    // fastMappedReadSeek must return an iterator, even if that iterator iterates through nothing.
                    return new NullSAMIterator(reads);
                } else {
                    readsTaken = readCount;
                    readsSeenAtLastPos = 0;
                    lastReadPos = GenomeLocParser.setStop(lastReadPos,-1);
                    CloseableIterator<SAMRecord> ret = iteratorPool.iterator(new MappedStreamSegment(lastReadPos));
                    return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, ret), readCount);
                }
            }
        }

        // if we're off the end of the last contig (into unmapped territory)
        if (rec != null && rec.getAlignmentStart() == 0) {
            readsTaken += readCount;
            intoUnmappedReads = true;
        }
        // else we're not off the end, store our location
        else if (rec != null) {
            int stopPos = rec.getAlignmentStart();
            if (stopPos < lastReadPos.getStart()) {
                lastReadPos = GenomeLocParser.createGenomeLoc(lastReadPos.getContigIndex() + 1, stopPos, stopPos);
            } else {
                lastReadPos = GenomeLocParser.setStart(lastReadPos,rec.getAlignmentStart());
            }
        }
        // in case we're run out of reads, get out
        else {
            throw new StingException("Danger: weve run out reads in fastMappedReadSeek");
            //return null;
        }
        bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);


        // return the iterator
        return bound;
    }

    /**
     * Even though the iterator has seeked to the correct location, there may be multiple reads at that location,
     * and we may have given some of them out already.  Move the iterator to the correct location using the readsAtLastPos variable
     *
     * @param iter the iterator
     */
    private void correctForReadPileupSeek( StingSAMIterator iter ) {
        // move the number of reads we read from the last pos
        boolean atLeastOneReadSeen = false; // we have a problem where some chomesomes don't have a single read (i.e. the chrN_random chrom.)
        while (iter.hasNext() && this.readsSeenAtLastPos > 0) {
            iter.next();
            --readsSeenAtLastPos;
            atLeastOneReadSeen = true;
        }
        if (readsSeenAtLastPos != 0 && atLeastOneReadSeen) {
            throw new SimpleDataSourceLoadException("Seek problem: reads at last position count != 0");
        }
    }


    /**
     * set the initial iterator
     *
     * @param readCount the number of reads
     * @param iter      the merging iterator
     *
     * @return a bounded read iterator at the first read position in the file.
     */
    private BoundedReadIterator InitialReadIterator( long readCount, CloseableIterator<SAMRecord> iter ) {
        BoundedReadIterator bound;
        bound = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads, iter), readCount);
        this.readsTaken = readCount;
        return bound;
    }

    /**
     * Filter reads based on user-specified criteria.
     *
     * @param enableVerification Verify the order of reads.
     * @param wrappedIterator the raw data source.
     * @param downsamplingFraction whether and how much to downsample the reads themselves (not at a locus).
     * @param filterZeroMappingQualityReads whether to filter zero mapping quality reads.
     * @param beSafeP Another trigger for the verifying iterator?  TODO: look into this.
     * @return An iterator wrapped with filters reflecting the passed-in parameters.  Will not be null.
     */
    private StingSAMIterator applyDecoratingIterators(boolean enableVerification,
                                                      StingSAMIterator wrappedIterator,
                                                      Double downsamplingFraction,
                                                      Boolean filterZeroMappingQualityReads,
                                                      Boolean beSafeP) {
        // NOTE: this (and other filtering) should be done before on-the-fly sorting
        //  as there is no reason to sort something that we will end of throwing away
        if (downsamplingFraction != null)
            wrappedIterator = new DownsampleIterator(wrappedIterator, downsamplingFraction);

        if (beSafeP != null && beSafeP && enableVerification)
            wrappedIterator = new VerifyingSamIterator(wrappedIterator);

        if ( filterZeroMappingQualityReads != null && filterZeroMappingQualityReads )
            wrappedIterator = StingSAMIteratorAdapter.adapt(wrappedIterator.getSourceInfo(),
                    new FilteringIterator(wrappedIterator, new ZeroMappingQualityReadFilterFunc()));

        return wrappedIterator;
    }

    private static class ZeroMappingQualityReadFilterFunc implements SamRecordFilter {
        public boolean filterOut(SAMRecord rec) {
            if (rec.getMappingQuality() == 0) {
                //System.out.printf("Filtering 0 mapping quality read %s%n", rec.format());
                return true;
            } else {
                return false;
            }
        }
    }
    
}

class SAMIteratorPool extends ResourcePool<ReadStreamPointer, StingSAMIterator> {
    /** Source information about the reads. */
    protected Reads reads;

    /**
     * A histogram of exactly what reads were removed from the input stream and why.
     */
    private SAMReadViolationHistogram violations = new SAMReadViolationHistogram();

    /** Is this a by-reads traversal or a by-locus? */
    protected boolean queryOverlapping;

    /** File header for the combined file. */
    protected SAMFileHeader header;

    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(SAMIteratorPool.class);

    public SAMIteratorPool( Reads reads ) {
        this.reads = reads;
        this.queryOverlapping = true;

        ReadStreamPointer streamPointer = createNewResource();
        this.header = streamPointer.getHeader();
        // Add this resource to the pool.
        this.addNewResource(streamPointer);
    }

    /** Get the combined header for all files in the iterator pool. */
    public SAMFileHeader getHeader() {
        return header;
    }

    /**
     * Returns a histogram of reads that were screened out, grouped by the nature of the error.
     * @return Histogram of reads.  Will not be null.
     */
    public SAMReadViolationHistogram getViolationHistogram() {
        return violations;
    }

    protected ReadStreamPointer selectBestExistingResource( DataStreamSegment segment, List<ReadStreamPointer> pointers ) {
        for (ReadStreamPointer pointer : pointers) {
            if (pointer.canAccessSegmentEfficiently(segment)) {
                return pointer;
            }
        }
        return null;
    }

    protected ReadStreamPointer createNewResource() {
        return new ReadStreamPointer(reads);
    }

    protected StingSAMIterator createIteratorFromResource( DataStreamSegment segment, ReadStreamPointer streamPointer ) {
        StingSAMIterator iterator = null;

        if (!queryOverlapping)
            iterator = streamPointer.getReadsContainedBy(segment);
        else {
            if (!( segment instanceof MappedStreamSegment ))
                throw new StingException("Segment is unmapped; true overlaps cannot be determined.");
            iterator = streamPointer.getReadsOverlapping((MappedStreamSegment) segment);
        }

        return new ReleasingIterator(new MalformedSAMFilteringIterator(header, iterator, violations));
    }

    protected void closeResource( ReadStreamPointer resource ) {
        resource.close();
    }

    private class ReleasingIterator implements StingSAMIterator {
        private final StingSAMIterator wrappedIterator;

        public Reads getSourceInfo() {
            return wrappedIterator.getSourceInfo();
        }

        public ReleasingIterator( StingSAMIterator wrapped ) {
            this.wrappedIterator = wrapped;
        }

        public ReleasingIterator iterator() {
            return this;
        }

        public void remove() {
            throw new UnsupportedOperationException("Can't remove from a StingSAMIterator");
        }

        public void close() {
            wrappedIterator.close();
            release(this);
        }

        public boolean hasNext() {
            return wrappedIterator.hasNext();
        }

        public SAMRecord next() {
            return wrappedIterator.next();
        }
    }

    public boolean isQueryOverlapping() {
        return queryOverlapping;
    }

    public void setQueryOverlapping( boolean queryOverlapping ) {
        this.queryOverlapping = queryOverlapping;
    }
}
