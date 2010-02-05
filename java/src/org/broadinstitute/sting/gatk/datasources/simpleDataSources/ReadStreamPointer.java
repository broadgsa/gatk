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

package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.SAMFileReader;

/**
 * Abstract class that models a current state in some category of reads.
 * @author hanna
 * @version 0.1
 */
abstract class ReadStreamPointer {
    /**
     * Describes the source of reads data.
     */
    protected final Reads sourceInfo;

    /**
     * Open handles to the reads info.
     */
    protected final SamFileHeaderMerger headerMerger;

    public ReadStreamPointer( Reads sourceInfo, SamFileHeaderMerger headerMerger ) {
        this.sourceInfo = sourceInfo;
        this.headerMerger = headerMerger;
    }

    /**
     * Can this pointer access the provided segment efficiently?
     * @param segment Segment to test.
     * @return True if it would be quick for this segment to access the given data.
     *         False if accessing this data would require some sort of reinitialization.
     */
    public abstract boolean canAccessSegmentEfficiently(DataStreamSegment segment);

    /**
     * Close this resource, destroying all file handles.
     */
    public void close() {
        for (SAMFileReader reader : headerMerger.getReaders())
            reader.close();
    }
    
    /**
     * Returns Reads data structure containing information about the reads data sources as well as
     * information about how they are downsampled, sorted, and filtered
     * @return
     */
    public Reads getReadsInfo() {
    	return sourceInfo;
    }

    /** 
     * Returns header merger: a class that keeps the mapping between original read groups and read groups
     * of the merged stream; merger also provides access to the individual file readers (and hence headers
     * too) maintained by the system. 
     * @return
     */
    public SamFileHeaderMerger getHeaderMerger() {
    	return headerMerger;
    }
    
    /**
     * Remove an iterator from service.
     * @param iterator The iterator to remove from service.  Must not be null.
     */
    public abstract void destroy( StingSAMIterator iterator );

    /**
     * Get a stream of all the reads that overlap a given segment.
     * @param segment Segment to check for overlaps.
     * @return An iterator over all reads overlapping the given segment.
     */
    public abstract StingSAMIterator getReadsOverlapping( DataStreamSegment segment );

    /**
     * Get a stream of all the reads that are completely contained by a given segment.
     * The segment can be mapped or unmapped.
     * @param segment Segment to check for containment..
     * @return An iterator over all reads contained by the given segment.
     */
    public abstract StingSAMIterator getReadsContainedBy( DataStreamSegment segment );
}

class MappedReadStreamPointer extends ReadStreamPointer {

    public MappedReadStreamPointer( Reads sourceInfo, SamFileHeaderMerger headerMerger ) {
        super( sourceInfo, headerMerger );
    }

    /**
     * MappedReadStreamPointers can access any segment efficiently.  Always return true.
     * @param segment Segment to test.
     * @return True.
     */
    public boolean canAccessSegmentEfficiently(DataStreamSegment segment) {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    public void destroy( StingSAMIterator iterator ) {
        iterator.close();
    }


    /**
     * {@inheritDoc}
     */
    @Override
    public StingSAMIterator getReadsOverlapping( DataStreamSegment segment ) {
        if(!(segment instanceof MappedStreamSegment))
            throw new UnsupportedOperationException("MappedReadStreamPointer cannot get reads overlapping an unmapped stream segment");
        MappedStreamSegment mappedSegment = (MappedStreamSegment)segment;

        MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );

        // The getStop() + 1 is a hack to work around an old bug in the way Picard created SAM files where queries
        // over a given interval would occasionally not pick up the last read in that interval.
        GenomeLoc bounds = mappedSegment.getBounds();
        mergingIterator.queryOverlapping( bounds.getContig(),
                                          (int)bounds.getStart(),
                                          (int)bounds.getStop()+ PlusOneFixIterator.PLUS_ONE_FIX_CONSTANT);

        return StingSAMIteratorAdapter.adapt(sourceInfo,mergingIterator);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public StingSAMIterator getReadsContainedBy( DataStreamSegment segment ) {
        if( !(segment instanceof MappedStreamSegment) )
            throw new StingException("Trying to access unmapped content from a mapped read stream pointer");
        MappedStreamSegment mappedSegment = (MappedStreamSegment)segment;
        MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );
        // NOTE: explicitly not using the queryOverlapping hack above since, according to the above criteria,
        //       we'd only miss reads that are one base long when performing a contained query.
        GenomeLoc bounds = mappedSegment.getBounds();
        mergingIterator.queryContained( bounds.getContig(),
                (int)bounds.getStart(),
                (int)bounds.getStop()+1);
        return StingSAMIteratorAdapter.adapt(sourceInfo,mergingIterator);
    }

    /**
     * Convert a mapped read stream pointer to an unmapped read stream pointer, transferring ownership
     * of the underlying file handles to the new container.
     * After doing this conversion, the source MappedReadStreamPointer should not be used.
     * @return
     */
    public UnmappedReadStreamPointer toUnmappedReadStreamPointer() {
        return new UnmappedReadStreamPointer( sourceInfo, headerMerger );
    }
}

class UnmappedReadStreamPointer extends ReadStreamPointer {
    /**
     * A pointer to the current position of this iterator in the read stream.
     */
    private PositionTrackingIterator unmappedIterator = null;

    public UnmappedReadStreamPointer( Reads sourceInfo, SamFileHeaderMerger headerMerger ) {
        super( sourceInfo, headerMerger );

        MergingSamRecordIterator2 mergingIterator = new MergingSamRecordIterator2( headerMerger, sourceInfo );
        mergingIterator.queryUnmappedReads();
        unmappedIterator = new PositionTrackingIterator( sourceInfo, mergingIterator, 0L );
    }

    /**
     * UnmappedReadStreamPointers are streams and can therefore access 'future' reads in the file quickly,
     * but reads already seen are difficult to seek to.
     * @param segment Segment to test.
     * @return True if this DataStreamSegment follows the current position.
     */
    public boolean canAccessSegmentEfficiently(DataStreamSegment segment) {
        if( !(segment instanceof UnmappedStreamSegment) )
            return false;
        return unmappedIterator.getPosition() <= ((UnmappedStreamSegment)segment).position;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public StingSAMIterator getReadsOverlapping( DataStreamSegment segment ) {
        throw new UnsupportedOperationException("Unable to determine overlapped reads of an unmapped segment");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public StingSAMIterator getReadsContainedBy( DataStreamSegment segment ) {
        if( !(segment instanceof UnmappedStreamSegment) )
            throw new StingException("Trying to access mapped content from an unmapped read stream pointer");

        UnmappedStreamSegment unmappedSegment = (UnmappedStreamSegment)segment;

        // Force the iterator to the next pending position.
        while(unmappedIterator.getPosition() < unmappedSegment.position)
            unmappedIterator.next();

        return new BoundedReadIterator(StingSAMIteratorAdapter.adapt(sourceInfo,unmappedIterator), unmappedSegment.size);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() {
        if( unmappedIterator != null )
            unmappedIterator.close();        
        super.close();
    }

    /**
     * {@inheritDoc}
     */
    public void destroy( StingSAMIterator iterator ) {
        // Don't destroy the iterator; reuse it.
    }

}

class EntireReadStreamPointer extends ReadStreamPointer {
    /**
     * Create a new pointer that can return info about the entire read stream.
     * @param sourceInfo Source info for the reads.
     * @param headerMerger Header merging apparatus.
     *
     */
    public EntireReadStreamPointer( Reads sourceInfo, SamFileHeaderMerger headerMerger ) {
        super( sourceInfo, headerMerger );
    }

    /**
     * An EntireReadStreamPointer can only efficiently access the entire file.
     * @param segment Segment to test.
     * @return true if requesting the entire stream.
     */
    public boolean canAccessSegmentEfficiently(DataStreamSegment segment) {
        return segment instanceof EntireStream;
    }

    /**
     * Get a stream of all the reads that overlap a given segment.
     * @param segment Segment to check for overlaps.
     * @return An iterator over all reads overlapping the given segment.
     */
    @Override
    public StingSAMIterator getReadsOverlapping( DataStreamSegment segment ) {
        if(!(segment instanceof EntireStream))
            throw new StingException("EntireReadStreamPointer can only get reads overlapping the entire stream.");
        return StingSAMIteratorAdapter.adapt(sourceInfo,new MergingSamRecordIterator2(headerMerger, sourceInfo));
    }

    /**
     * Get a stream of all the reads that are completely contained by a given segment.
     * @param segment Segment to check for containment..
     * @return An iterator over all reads contained by the given segment.
     */
    @Override
    public StingSAMIterator getReadsContainedBy( DataStreamSegment segment ) {
        if(!(segment instanceof EntireStream))
            throw new StingException("EntireReadStreamPointer can only get reads contained by the entire stream.");
        return StingSAMIteratorAdapter.adapt(sourceInfo,new MergingSamRecordIterator2(headerMerger, sourceInfo));        
    }

    /**
     * {@inheritDoc}
     */
    public void destroy( StingSAMIterator iterator ) {
        iterator.close();
    }    

}