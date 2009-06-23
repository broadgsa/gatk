package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.utils.StingException;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.util.List;
import java.util.ArrayList;
import java.io.File;
/**
 * User: hanna
 * Date: Jun 23, 2009
 * Time: 6:49:04 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Maintains a pointer into a stream of reads.  Tracks state between mapped and unmapped.
 * For mapped, assumes that the user will query directly to where they want; closes the iterator after each use.
 * For unmapped, assumes that the user will walk through the entire stream.  Keeps the iterator open permanently.
 */
enum MappingType { MAPPED, UNMAPPED }

class ReadStreamPointer {
    /** our log, which we want to capture anything from this class */
    protected static Logger logger = Logger.getLogger(ReadStreamPointer.class);

    /**
     * Describes the source of reads data.
     */
    private final Reads sourceInfo;

    /**
     * Open handles to the reads info.
     */
    private final SamFileHeaderMerger headerMerger;

    /**
     * The (possibly merged) header for the input fileset.
     */
    private final SAMFileHeader header;

    /**
     * In which bucket of reads does this pointer live?
     */
    private MappingType streamPosition = MappingType.MAPPED;

    /**
     * A pointer to the current position of this iterator in the read stream.
     */
    private PositionTrackingIterator unmappedIterator = null;

    public ReadStreamPointer( Reads sourceInfo ) {
        this.sourceInfo = sourceInfo;
        this.headerMerger = createHeaderMerger(sourceInfo, SAMFileHeader.SortOrder.coordinate);
        this.header = this.headerMerger.getMergedHeader();
    }

    /**
     * Gets the header information for the read stream.
     * @return Header information for the read stream.
     */
    public SAMFileHeader getHeader() {
        return header;
    }

    /**
     * Can this pointer be efficiently used to access the given segment?
     * @param segment Segment to inspect.
     * @return True if the segment can be accessed efficiently, false otherwise.
     */
    public boolean canAccessSegmentEfficiently( DataStreamSegment segment ) {
        switch( streamPosition ) {
            case MAPPED:
                return true;
            case UNMAPPED:
                if( segment instanceof MappedStreamSegment )
                    return false;
                else if( segment instanceof UnmappedStreamSegment ) {
                    UnmappedStreamSegment unmappedSegment = (UnmappedStreamSegment)segment;
                    return unmappedIterator.position <= unmappedSegment.position;
                }
                else
                    throw new StingException("Unsupported stream segment type: " + segment.getClass());
            default:
                throw new StingException("Pointer has hit illegal stream position; current position is " + streamPosition);

        }
    }

    public void close() {
        if( unmappedIterator != null )
            unmappedIterator.close();
        for (SAMFileReader reader : headerMerger.getReaders())
            reader.close();
    }

    /**
     * Get a stream of all the reads that overlap a given segment.
     * @param segment Segment to check for overlaps.
     * @return An iterator over all reads overlapping the given segment.
     */
    public StingSAMIterator getReadsOverlapping( MappedStreamSegment segment ) {
        MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );
        mergingIterator.queryOverlapping( segment.locus.getContig(),
                                          (int)segment.locus.getStart(),
                                          (int)segment.locus.getStop());
        return StingSAMIteratorAdapter.adapt(sourceInfo,mergingIterator);
    }

    public StingSAMIterator getReadsContainedBy( DataStreamSegment segment ) {
        if( segment instanceof MappedStreamSegment ) {
            MappedStreamSegment mappedSegment = (MappedStreamSegment)segment;
            MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );
            mergingIterator.queryContained( mappedSegment.locus.getContig(),
                                            (int)mappedSegment.locus.getStart(),
                                            (int)mappedSegment.locus.getStop());
            return StingSAMIteratorAdapter.adapt(sourceInfo,mergingIterator);
        }
        else if( segment instanceof UnmappedStreamSegment ) {
            UnmappedStreamSegment unmappedSegment = (UnmappedStreamSegment)segment;

            // If the stream position has not flipped over to the unmapped state, do some initialization.
            if( streamPosition == MappingType.MAPPED ) {
                MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );
                mergingIterator.queryUnmappedReads();
                unmappedIterator = new PositionTrackingIterator( sourceInfo, mergingIterator, 0L );
                streamPosition = MappingType.UNMAPPED;
            }
            else {
                if( streamPosition != MappingType.UNMAPPED || unmappedIterator == null )
                    throw new StingException("Illegal state: iterator has fetched all mapped reads but has not properly transition to unmapped reads");

                // Force the iterator to the next pending position.
                while(unmappedIterator.position < unmappedSegment.position)
                    unmappedIterator.next();
            }

            return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(sourceInfo,unmappedIterator), unmappedSegment.size);
        }
        else
            throw new StingException("Unable to handle stream segment of type" + segment.getClass());
    }

    /**
     * A private function that, given the internal file list, generates a merging construct for
     * all available files.
     * @param reads source information about the reads.
     * @param SORT_ORDER sort order for the reads.
     * @return a list of SAMFileReaders that represent the stored file names
     */
    protected SamFileHeaderMerger createHeaderMerger( Reads reads, SAMFileHeader.SortOrder SORT_ORDER )
            throws SimpleDataSourceLoadException {
        // right now this is pretty damn heavy, it copies the file list into a reader list every time
        List<SAMFileReader> lst = new ArrayList<SAMFileReader>();
        for (File f : reads.getReadsFiles()) {
            SAMFileReader reader = new SAMFileReader(f, true);
            reader.setValidationStringency(reads.getValidationStringency());

            final SAMFileHeader header = reader.getFileHeader();
            logger.debug(String.format("Sort order is: " + header.getSortOrder()));

            if (reader.getFileHeader().getReadGroups().size() < 1) {
                //logger.warn("Setting header in reader " + f.getName());
                SAMReadGroupRecord rec = new SAMReadGroupRecord(f.getName());
                rec.setLibrary(f.getName());
                rec.setSample(f.getName());

                reader.getFileHeader().addReadGroup(rec);
            }

            lst.add(reader);
        }
        return new SamFileHeaderMerger(lst,SORT_ORDER,true);
    }

    private class PositionTrackingIterator implements StingSAMIterator {
        /**
         * Source information about the reads.
         */
        private Reads sourceInfo;

        /**
         * The iterator being tracked.
         */
        private CloseableIterator<SAMRecord> iterator;

        /**
         * Current position within the tracked iterator.
         */
        private long position;

        /**
         * {@inheritDoc}
         */
        public Reads getSourceInfo() {
            return sourceInfo;
        }

        /**
         * Retrieves the current position of the iterator.  The 'current position' of the iterator is defined as
         * the coordinate of the read that will be returned if next() is called.
         * @return The current position of the iterator.
         */
        public long getPosition() {
            return position;
        }

        /**
         * Create a new iterator wrapping the given position, assuming that the reader is <code>position</code> reads
         * into the sequence.
         * @param sourceInfo Information about where these reads came from.
         * @param iterator Iterator to wraps.
         * @param position Non-negative position where the iterator currently sits.
         */
        public PositionTrackingIterator( Reads sourceInfo, CloseableIterator<SAMRecord> iterator, long position ) {
            this.sourceInfo = sourceInfo;
            this.iterator = iterator;
            this.position = position;
        }

        /**
         * {@inheritDoc}
         */
        public boolean hasNext() {
            return iterator.hasNext();
        }

        /**
         * Try to get the next read in the list.  If a next read is available, increment the position.
         * @return next read in the list, if available.
         */
        public SAMRecord next() {
            try {
                return iterator.next();
            }
            finally {
                position++;
            }
        }

        /**
         * {@inheritDoc}
         */
        public StingSAMIterator iterator() {
            return this;
        }

        /**
         * {@inheritDoc}
         */
        public void close() {
            // Position tracking iterators are constant through the life of the traversal.  Don't close them.
            // TODO: This is an artifact of the fact that pooled query iterators need to be closed, but pooled unmapped
            // TODO: iterators must not be.  Clean this up!
        }

        /**
         * {@inheritDoc}
         */
        public void remove() { throw new UnsupportedOperationException("Cannot remove from a StingSAMIterator"); }

    }
}
