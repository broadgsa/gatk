package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.sam.SamFileHeaderMerger;

import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.MonolithicShard;
import org.broadinstitute.sting.gatk.datasources.shards.ReadDelimitedReadShard;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.SAMReadViolationHistogram;

import java.util.*;
import java.io.File;

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
public class IndexDrivenSAMDataSource extends SAMDataSource {
    // used for the reads case, the last count of reads retrieved
    long readsTaken = 0;

    // our last genome loc position
    protected GenomeLoc lastReadPos = null;

    // do we take unmapped reads
    private boolean includeUnmappedReads = true;

    // reads based traversal variables
    private boolean intoUnmappedReads = false;
    private int readsSeenAtLastPos = 0;

    /**
     * A histogram of exactly what reads were removed from the input stream and why.
     */
    private SAMReadViolationHistogram violations = new SAMReadViolationHistogram();

    // A pool of SAM iterators.
    private SAMResourcePool resourcePool = null;

    private GenomeLoc mLastInterval = null;

    /**
     * Returns a histogram of reads that were screened out, grouped by the nature of the error.
     * @return Histogram of reads.  Will not be null.
     */
    public SAMReadViolationHistogram getViolationHistogram() {
        return violations;
    }

    /**
     * constructor, given sam files
     *
     * @param reads the list of sam files
     */
    public IndexDrivenSAMDataSource( Reads reads ) throws SimpleDataSourceLoadException {
        super(reads);
        resourcePool = new SAMResourcePool(reads);
    }

    /**
     * Do all BAM files backing this data source have an index?  The case where hasIndex() is false
     * is supported, but only in a few extreme cases.
     * @return True if an index is present; false otherwise.
     */
    public boolean hasIndex() {
        return resourcePool.hasIndex;
    }

    /**
     * Gets the (potentially merged) SAM file header.
     *
     * @return SAM file header.
     */
    public SAMFileHeader getHeader() {
        return resourcePool.getHeader();
    }

    /**
     * Returns a mapping from original input files to the SAMFileReaders
     *
     * @return the mapping
     */
    public Map<File, SAMFileReader> getFileToReaderMapping() {
        return resourcePool.getFileToReaderMapping();
    }

    /**
     * Returns Reads data structure containing information about the reads data sources placed in this pool as well as
     * information about how they are downsampled, sorted, and filtered
     * @return
     */
    public Reads getReadsInfo() { return reads; }

    /**
     * Returns header merger: a class that keeps the mapping between original read groups and read groups
     * of the merged stream; merger also provides access to the individual file readers (and hence headers
     * prior to the merging too) maintained by the system.
     * @return
     */
    public Collection<SAMFileReader> getReaders() { return resourcePool.getHeaderMerger().getReaders(); }

    /** Returns true if there are read group duplicates within the merged headers. */
    public boolean hasReadGroupCollisions() {
        return resourcePool.getHeaderMerger().hasReadGroupCollisions();
    }

    /** Returns the read group id that should be used for the input read and RG id. */
    public String getReadGroupId(final SAMFileReader reader, final String originalReadGroupId) {
        return resourcePool.getHeaderMerger().getReadGroupId(reader,originalReadGroupId);
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
        resourcePool.setQueryOverlapping(queryOverlapping);

        StingSAMIterator iterator = null;
        if (shard.getShardType() == Shard.ShardType.READ) {
            iterator = seekRead(shard);
            iterator = applyDecoratingIterators(true,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                    reads.getSupplementalFilters());
        } else if (shard.getShardType() == Shard.ShardType.LOCUS) {
            iterator = seekLocus(shard);
            iterator = applyDecoratingIterators(false,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                    reads.getSupplementalFilters());
        } else if ((shard.getShardType() == Shard.ShardType.LOCUS_INTERVAL) ||
                   (shard.getShardType() == Shard.ShardType.READ_INTERVAL)) {
            iterator = seekLocus(shard);
            iterator = applyDecoratingIterators(false,
                    iterator,
                    reads.getDownsamplingFraction(),
                    reads.getValidationExclusionList().contains(ValidationExclusion.TYPE.NO_READ_ORDER_VERIFICATION),
                    reads.getSupplementalFilters());

            // add the new overlapping detection iterator, if we have a last interval and we're a read based shard
            if(shard.getGenomeLocs().size() > 1)
                throw new StingException("This SAMDataSource does not support multiple intervals within a single shard");
            GenomeLoc shardGenomeLoc = shard.getGenomeLocs().get(0);
            if (mLastInterval != null && shard.getShardType() == Shard.ShardType.READ_INTERVAL )
                iterator = new PlusOneFixIterator(shardGenomeLoc,new IntervalOverlapIterator(iterator,mLastInterval,false));
            mLastInterval = shardGenomeLoc;
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
     * @param shard the shard containing the genome location to extract data for
     *
     * @return an iterator for that region
     */
    private StingSAMIterator seekLocus( Shard shard ) throws SimpleDataSourceLoadException {
        if(shard instanceof MonolithicShard)
            return createIterator(new EntireStream());

        if( getHeader().getSequenceDictionary().getSequences().size() == 0 )
            throw new StingException("Unable to seek to the given locus; reads data source has no alignment information.");

        if(shard.getGenomeLocs().size() > 1)
            throw new StingException("This SAMDataSource does not support multiple intervals within a single shard");
        GenomeLoc shardGenomeLoc = shard.getGenomeLocs().get(0);

        return createIterator( new MappedStreamSegment(Collections.singletonList(shardGenomeLoc)) );
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
    private StingSAMIterator seekRead( Shard shard ) throws SimpleDataSourceLoadException {
        if(shard instanceof MonolithicShard)
            return createIterator(new EntireStream());

        ReadDelimitedReadShard readShard = (ReadDelimitedReadShard)shard;
        StingSAMIterator iter = null;

        // If there are no entries in the sequence dictionary, there can't possibly be any unmapped reads.  Force state to 'unmapped'.
        if( isSequenceDictionaryEmpty() )
            intoUnmappedReads = true;

        if (!intoUnmappedReads) {
            if (lastReadPos == null) {
                lastReadPos = GenomeLocParser.createGenomeLoc(getHeader().getSequenceDictionary().getSequence(0).getSequenceIndex(), 0, Integer.MAX_VALUE);
                iter = createIterator(new MappedStreamSegment(Collections.singletonList(lastReadPos)));
                return InitialReadIterator(readShard.getSize(), iter);
            } else {
                lastReadPos = GenomeLocParser.setStop(lastReadPos,-1);
                iter = fastMappedReadSeek(readShard.getSize(), StingSAMIteratorAdapter.adapt(reads, createIterator(new MappedStreamSegment(Collections.singletonList(lastReadPos)))));
            }

            if (intoUnmappedReads && !includeUnmappedReads)
                readShard.signalDone();
        }

        if (intoUnmappedReads && includeUnmappedReads) {
            if (iter != null)
                iter.close();
            iter = toUnmappedReads(readShard.getSize());
            if (!iter.hasNext())
                readShard.signalDone();
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
     * @param resourcePool Custom mock iterator pool.
     */
    void setResourcePool( SAMResourcePool resourcePool ) {
        this.resourcePool = resourcePool;
    }

    /**
     * Retrieve unmapped reads.
     *
     * @param readCount how many reads to retrieve
     *
     * @return the bounded iterator that you can use to get the intervaled reads from
     */
    StingSAMIterator toUnmappedReads( long readCount ) {
        StingSAMIterator iter = createIterator(new UnmappedStreamSegment(readsTaken, readCount));
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

        // Assuming that lastReadPos should never be null, because this is a mappedReadSeek
        // and initial queries are handled by the previous conditional.
        int lastContig = lastReadPos.getContigIndex();
        int lastPos = (int)lastReadPos.getStart();

        while (x < readsTaken) {
            if (iter.hasNext()) {
                rec = iter.next();
                if (lastContig == rec.getReferenceIndex() && lastPos == rec.getAlignmentStart()) ++this.readsSeenAtLastPos;
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
                    CloseableIterator<SAMRecord> ret = createIterator(new MappedStreamSegment(Collections.singletonList(lastReadPos)));
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
     * Determines whether the BAM file is completely unsequenced.  Requires that the resource pool be initialized.
     * @return True if the sequence dictionary is completely empty.  False otherwise.
     */
    private boolean isSequenceDictionaryEmpty() {
        return getHeader().getSequenceDictionary().isEmpty();
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
        for(int i = 0; i < this.readsSeenAtLastPos && iter.hasNext(); i++,iter.next())
            atLeastOneReadSeen = true;
        if (readsSeenAtLastPos > 0 && !atLeastOneReadSeen) {
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
     * Creates an iterator over the selected segment, from a resource pulled from the pool.
     * @param segment Segment over which to gather reads.
     * @return An iterator over just the reads in the given segment.
     */
    private StingSAMIterator createIterator( DataStreamSegment segment ) {
        StingSAMIterator iterator = resourcePool.iterator(segment);
        StingSAMIterator malformedWrappedIterator =  new MalformedSAMFilteringIterator( getHeader(), iterator, violations );
        StingSAMIterator readWrappingIterator = new ReadWrappingIterator(malformedWrappedIterator);
        return readWrappingIterator;
    }
}


