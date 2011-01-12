package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.FlashBackIterator;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.reflect.Type;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 10:04:12 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A data source which provides a single type of reference-ordered data.
 */
public class ReferenceOrderedDataSource implements SimpleDataSource {
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
        return builder.getAvailableTrackNamesAndTypes().get(fileDescriptor.getType().toUpperCase());                
    }

    public Class getRecordType() {
        return builder.createCodec(getType(),getName()).getFeatureType();
    }

    public File getFile() {
        return new File(fileDescriptor.getFile());
    }

    public Object getHeader() {
        return header;    
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
     * @param shard Shard that points to the selected position.
     * @return Iterator through the data.
     */
    public LocationAwareSeekableRODIterator seek( Shard shard ) {
        DataStreamSegment dataStreamSegment = shard.getGenomeLocs().size() != 0 ? new MappedStreamSegment(shard.getGenomeLocs().get(0)) : new EntireStream();
        return iteratorPool.iterator(dataStreamSegment);
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
    public ReferenceOrderedDataPool(RMDTriplet fileDescriptor,RMDTrackBuilder builder,SAMSequenceDictionary sequenceDictionary,GenomeLocParser genomeLocParser,boolean flashbackData) {
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
            throw new ReviewedStingException("BUG: Tried to create multiple iterators over streaming ROD interface");
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
            throw new ReviewedStingException("Unable to find a ROD iterator for segments of type " + segment.getClass());
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
        } catch (IOException e) {
            throw new ReviewedStingException("Unable to create iterator for rod named " + fileDescriptor.getName(),e);
        }
    }

    @Override
    protected void closeResource(RMDTrack track) {
        track.close();
    }
}


