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

package org.broadinstitute.gatk.engine.datasources.rmd;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.gatk.utils.refdata.SeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrack;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.FlashBackIterator;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.List;

/**
 * A pool of reference-ordered data iterators.
 */
class ReferenceOrderedDataPool extends ResourcePool<LocationAwareSeekableRODIterator, LocationAwareSeekableRODIterator> {
    // the reference-ordered data itself.
    private final RMDTriplet fileDescriptor;

    // our tribble track builder
    private final RMDTrackBuilder builder;

    /**
     * The header from this RMD, if present.
     */
    private final Object header;

    /**
     * The sequence dictionary from this ROD.  If no sequence dictionary is present, this dictionary will be the same as the reference's.
     */
    private final SAMSequenceDictionary sequenceDictionary;

    boolean flashbackData = false;
    public ReferenceOrderedDataPool(RMDTriplet fileDescriptor,RMDTrackBuilder builder,SAMSequenceDictionary sequenceDictionary, GenomeLocParser genomeLocParser,boolean flashbackData) {
        super(sequenceDictionary,genomeLocParser);
        this.fileDescriptor = fileDescriptor;
        this.builder = builder;
        this.flashbackData = flashbackData;

        // prepopulate one RMDTrack
        LocationAwareSeekableRODIterator iterator = createNewResource();
        this.addNewResource(iterator);

        // Pull the proper header and sequence dictionary from the prepopulated track.
        this.header = iterator.getHeader();
        this.sequenceDictionary = iterator.getSequenceDictionary();
    }

    /**
     * Gets the header used by this resource pool.
     * @return Header used by this resource pool.
     */
    public Object getHeader() {
        return header;
    }

    /**
     * Gets the sequence dictionary built into the ROD index file.
     * @return Sequence dictionary from the index file.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    /**
     * Create a new iterator from the existing reference-ordered data.  This new iterator is expected
     * to be completely independent of any other iterator.
     * @return The newly created resource.
     */
    public LocationAwareSeekableRODIterator createNewResource() {
        if(numIterators() > 0)
            throw new ReviewedGATKException("BUG: Tried to create multiple iterators over streaming ROD interface");
        RMDTrack track = builder.createInstanceOfTrack(fileDescriptor);
        LocationAwareSeekableRODIterator iter = new SeekableRODIterator(track.getHeader(),track.getSequenceDictionary(),referenceSequenceDictionary,genomeLocParser,track.getIterator());
        return (flashbackData) ? new FlashBackIterator(iter) : iter;
    }

    /**
     * Finds the best existing ROD iterator from the pool.  In this case, the best existing ROD is defined as
     * the first one encountered that is at or before the given position.
     * @param segment @{inheritedDoc}
     * @param resources @{inheritedDoc}
     * @return @{inheritedDoc}
     */
    public LocationAwareSeekableRODIterator selectBestExistingResource( DataStreamSegment segment, List<LocationAwareSeekableRODIterator> resources ) {
        if(segment instanceof MappedStreamSegment) {
            GenomeLoc position = ((MappedStreamSegment)segment).getLocation();

            for( LocationAwareSeekableRODIterator RODIterator : resources ) {

                if( (RODIterator.position() == null && RODIterator.hasNext()) ||
                    (RODIterator.position() != null && RODIterator.position().isBefore(position)) )
                    return RODIterator;
                if (RODIterator.position() != null && RODIterator instanceof FlashBackIterator && ((FlashBackIterator)RODIterator).canFlashBackTo(position)) {
                    ((FlashBackIterator)RODIterator).flashBackTo(position);
                    return RODIterator;
                }

            }
            return null;
        }
        else if(segment instanceof EntireStream) {
            // Asking for a segment over the entire stream, so by definition, there is no best existing resource.
            // Force the system to create a new one.
            return null;
        }
        else {
            throw new ReviewedGATKException("Unable to find a ROD iterator for segments of type " + segment.getClass());
        }
    }

    /**
     * In this case, the iterator is the resource.  Pass it through.
     */
    public LocationAwareSeekableRODIterator createIteratorFromResource( DataStreamSegment segment, LocationAwareSeekableRODIterator resource ) {
        return resource;
    }

    /**
     * kill the buffers in the iterator
     */
    public void closeResource( LocationAwareSeekableRODIterator resource ) {
        if (resource instanceof FlashBackIterator) ((FlashBackIterator)resource).close();
    }
}
