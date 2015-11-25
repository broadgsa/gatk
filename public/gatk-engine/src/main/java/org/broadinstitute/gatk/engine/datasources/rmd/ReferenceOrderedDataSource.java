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
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.refdata.SeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrack;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RMDTriplet;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.Type;
import java.util.List;

/**
 * A data source which provides a single type of reference-ordered data.
 */
public class ReferenceOrderedDataSource {
    /**
     * The reference-ordered data itself.
     */
    private final RMDTriplet fileDescriptor;

    /**
     * The header associated with this VCF, if any.
     */
    private final Object header;

    /**
     * The private sequence dictionary associated with this RMD.
     */
    private final SAMSequenceDictionary sequenceDictionary;

    /**
     * The builder to use when constructing new reference-ordered data readers.
     */
    private final RMDTrackBuilder builder;

    /**
     * A pool of iterators for navigating through the genome.
     */
    private final ResourcePool<?,LocationAwareSeekableRODIterator> iteratorPool;

    /**
     * Create a new reference-ordered data source.
     */
    public ReferenceOrderedDataSource(RMDTriplet fileDescriptor,
                                      RMDTrackBuilder builder,
                                      SAMSequenceDictionary referenceSequenceDictionary,
                                      GenomeLocParser genomeLocParser,
                                      boolean flashbackData ) {
        this.fileDescriptor = fileDescriptor;
        this.builder = builder;

        // TODO: Unify the two blocks of code below by creating a ReferenceOrderedDataPool base class of a coherent type (not RMDTrack for one and SeekableIterator for the other).
        if (fileDescriptor.getStorageType() != RMDTriplet.RMDStorageType.STREAM) {
            iteratorPool = new ReferenceOrderedQueryDataPool(fileDescriptor,
                                                             builder,
                                                             referenceSequenceDictionary,
                                                             genomeLocParser);
            this.header = ((ReferenceOrderedQueryDataPool)iteratorPool).getHeader();
            this.sequenceDictionary = ((ReferenceOrderedQueryDataPool)iteratorPool).getSequenceDictionary();
        }
        else {
            iteratorPool = new ReferenceOrderedDataPool(fileDescriptor,
                                                        builder,
                                                        referenceSequenceDictionary,
                                                        genomeLocParser,
                                                        flashbackData);
            this.header = ((ReferenceOrderedDataPool)iteratorPool).getHeader();
            this.sequenceDictionary = ((ReferenceOrderedDataPool)iteratorPool).getSequenceDictionary();
        }
    }

    /**
     * Return the name of the underlying reference-ordered data.
     * @return Name of the underlying rod.
     */
    public String getName() {
        return fileDescriptor.getName();
    }

    public Class getType() {
        return builder.getFeatureManager().getByTriplet(fileDescriptor).getCodecClass();
    }

    public Class getRecordType() {
        return builder.getFeatureManager().getByTriplet(fileDescriptor).getFeatureClass();
    }

    public File getFile() {
        return new File(fileDescriptor.getFile());
    }

    public Object getHeader() {
        return header;
    }

    public Tags getTags() {
        return fileDescriptor.getTags();
    }
    
    public String getTagValue( final String key ) {
        return fileDescriptor.getTags().getValue( key );
    }


    /**
     * Retrieves the sequence dictionary created by this ROD.
     * @return
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    /**
     * helper function for determining if we are the same track based on name and record type
     *
     * @param name the name to match
     * @param type the type to match
     *
     * @return true on a match, false if the name or type is different
     */
    public boolean matchesNameAndRecordType(String name, Type type) {
        return (name.equals(fileDescriptor.getName()) && (type.getClass().isAssignableFrom(getType().getClass())));
    }

    /**
     * Seek to the specified position and return an iterator through the data.
     *
     * @param loc GenomeLoc that points to the selected position.
     *
     * @return Iterator through the data.
     */
    public LocationAwareSeekableRODIterator seek(GenomeLoc loc) {
        DataStreamSegment dataStreamSegment = loc != null ? new MappedStreamSegment(loc) : new EntireStream();
        return iteratorPool.iterator(dataStreamSegment);
    }


    /**
     * Close the specified iterator, returning it to the pool.
     * @param iterator Iterator to close.
     */
    public void close( LocationAwareSeekableRODIterator iterator ) {
        iteratorPool.release(iterator);
    }

}

/**
 * a data pool for the new query based RODs
 */
class ReferenceOrderedQueryDataPool extends ResourcePool<RMDTrack,LocationAwareSeekableRODIterator> {
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

    public ReferenceOrderedQueryDataPool(RMDTriplet fileDescriptor, RMDTrackBuilder builder, SAMSequenceDictionary referenceSequenceDictionary, GenomeLocParser genomeLocParser) {
        super(referenceSequenceDictionary,genomeLocParser);
        this.fileDescriptor = fileDescriptor;
        this.builder = builder;

        // prepopulate one RMDTrack
        RMDTrack track = builder.createInstanceOfTrack(fileDescriptor);
        this.addNewResource(track);

        // Pull the proper header and sequence dictionary from the prepopulated track.
        this.header = track.getHeader();
        this.sequenceDictionary = track.getSequenceDictionary();
    }

    public Object getHeader() {
        return header;
    }

    public SAMSequenceDictionary getSequenceDictionary() {
        return sequenceDictionary;
    }

    @Override
    protected RMDTrack createNewResource() {
        return builder.createInstanceOfTrack(fileDescriptor);
    }

    @Override
    protected RMDTrack selectBestExistingResource(DataStreamSegment segment, List<RMDTrack> availableResources) {
        for (RMDTrack reader : availableResources)
            if (reader != null) return reader;
        return null;
    }

    @Override
    protected LocationAwareSeekableRODIterator createIteratorFromResource(DataStreamSegment position, RMDTrack track) {
        try {
            if (position instanceof MappedStreamSegment) {
                GenomeLoc pos = ((MappedStreamSegment) position).locus;
                return new SeekableRODIterator(header,sequenceDictionary,referenceSequenceDictionary,genomeLocParser,track.query(pos));
            } else {
                return new SeekableRODIterator(header,sequenceDictionary,referenceSequenceDictionary,genomeLocParser,track.getIterator());
            }
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotReadInputFile(fileDescriptor.getName(), "it could not be found");
        } catch (IOException | RuntimeException e) {
            throw new ReviewedGATKException("Unable to create iterator for rod named " + fileDescriptor.getName(),e);
        }

    }

    @Override
    protected void closeResource(RMDTrack track) {
        track.close();
    }
}