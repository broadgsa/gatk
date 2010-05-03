package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.QueryableTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.utils.FlashBackIterator;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
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
    private final RMDTrack rod;

    /**
     * A pool of iterators for navigating through the genome.
     */
    private final ReferenceOrderedDataPool iteratorPool;

    /**
     * Create a new reference-ordered data source.
     * @param rod the reference ordered data
     */
    public ReferenceOrderedDataSource( Walker walker, RMDTrack rod) {
        this.rod = rod;
        if (rod.supportsQuery()) iteratorPool = null;
        else iteratorPool = new ReferenceOrderedDataPool( walker, rod );
    }

    /**
     * Return the name of the underlying reference-ordered data.
     * @return Name of the underlying rod.
     */
    public String getName() {
        return this.rod.getName(); 
    }

    /**
     * Return the underlying reference-ordered data.
     * @return the underlying rod.
     */
    public RMDTrack getReferenceOrderedData() {
        return this.rod;
    }

    /**
     * Seek to the specified position and return an iterator through the data.
     * @param shard Shard that points to the selected position.
     * @return Iterator through the data.
     */
    public LocationAwareSeekableRODIterator seek( Shard shard ) {
        if (iteratorPool == null) // use query
            return getQuery(shard.getGenomeLocs() == null || shard.getGenomeLocs().size() == 0 ? null : shard.getGenomeLocs());
        DataStreamSegment dataStreamSegment = shard.getGenomeLocs().size() != 0 ? new MappedStreamSegment(shard.getGenomeLocs().get(0)) : new EntireStream();
        LocationAwareSeekableRODIterator RODIterator = iteratorPool.iterator(dataStreamSegment);
        return RODIterator;
    }

    /**
     * Seek to the specified position and return an iterator through the data.
     *
     * @param loc GenomeLoc that points to the selected position.
     *
     * @return Iterator through the data.
     */
    public LocationAwareSeekableRODIterator seek(GenomeLoc loc) {
        if (iteratorPool == null) // use query
            return getQuery(loc == null ? null : Arrays.asList(loc));
        DataStreamSegment dataStreamSegment = loc != null ? new MappedStreamSegment(loc) : new EntireStream();
        LocationAwareSeekableRODIterator RODIterator = iteratorPool.iterator(dataStreamSegment);
        return RODIterator;
    }

    /**
     * assuming the ROD is a queryable ROD, use that interface to get an iterator to the selected region
     * @param loc the region to query for
     * @return a LocationAwareSeekableRODIterator over the selected region
     */
    private LocationAwareSeekableRODIterator getQuery(List<GenomeLoc> loc) {
        if (loc == null) // for the mono shard case
            return new SeekableRODIterator(rod.getIterator());
        return new StitchingLocationAwareSeekableRODIterator(loc,(QueryableTrack)rod);
    }

    /**
     * Close the specified iterator, returning it to the pool.
     * @param iterator Iterator to close.
     */
    public void close( LocationAwareSeekableRODIterator iterator ) {
        if (iteratorPool != null) iteratorPool.release(iterator);
    }

}

/**
 * A pool of reference-ordered data iterators.
 */
class ReferenceOrderedDataPool extends ResourcePool<LocationAwareSeekableRODIterator, LocationAwareSeekableRODIterator> {
    private final RMDTrack rod;
    boolean flashbackData = false;
    public ReferenceOrderedDataPool( Walker walker, RMDTrack rod ) {
        if (walker instanceof ReadWalker) flashbackData = true; // && (rod.getType() != IntervalRod.class)
        this.rod = rod;
    }

    /**
     * Create a new iterator from the existing reference-ordered data.  This new iterator is expected
     * to be completely independent of any other iterator.
     * @return The newly created resource.
     */
    public LocationAwareSeekableRODIterator createNewResource() {
        LocationAwareSeekableRODIterator iter = new SeekableRODIterator(rod.getIterator());
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
            GenomeLoc position = ((MappedStreamSegment)segment).getFirstLocation();

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
            throw new StingException("Unable to find a ROD iterator for segments of type " + segment.getClass());
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
 * stitch together the multiple calls to seek (since shards can have multiple intervals now)
 * on the underlying Tribble track into one seamless iteration
 */
class StitchingLocationAwareSeekableRODIterator implements LocationAwareSeekableRODIterator {

    // the list of intervals we're iterating over
    private final LinkedList<GenomeLoc> locationList;

    // The reference-ordered data itself.
    private final QueryableTrack rod;

    // the current iterator
    private SeekableRODIterator iterator;

    StitchingLocationAwareSeekableRODIterator(List<GenomeLoc> list, QueryableTrack rmd) {
        rod = rmd;
        locationList = new LinkedList<GenomeLoc>();
        locationList.addAll(list);
        fetchNextInterval();
    }

    @Override
    public GenomeLoc peekNextLocation() {
        if (iterator == null) return null;
        return iterator.peekNextLocation();
    }

    @Override
    public GenomeLoc position() {
        if (iterator == null) return null;
        return iterator.position();
    }

    @Override
    public RODRecordList seekForward(GenomeLoc interval) {
        RODRecordList list = iterator.seekForward(interval);
        if (list == null) { // we were unable to seek the current interval to the location
            fetchNextInterval();
            list = iterator.seekForward(interval);
        }
        return list;
    }

    @Override
    public boolean hasNext() {
        if (iterator == null) return false;
        return iterator.hasNext();
    }

    @Override
    public RODRecordList next() {
        if (!hasNext()) throw new IllegalStateException("StitchingLocationAwareSeekableRODIterator: We do not have a next");
        RODRecordList list = iterator.next();
        if (!iterator.hasNext()) fetchNextInterval();
        return list;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("\"Thou shall not remove()!\" - Software Engineering Team");
    }

    private void fetchNextInterval() {
        if (locationList != null && locationList.size() > 0) {
            GenomeLoc loc = locationList.getFirst();
            locationList.removeFirst();
            if (rod == null) throw new StingException("Unable to query(), target rod is null, next location = " + ((locationList != null) ? locationList.getFirst() : "null"));
            try {
                iterator = new SeekableRODIterator(rod.query(loc));
            } catch (IOException e) {
                throw new StingException("Unable to query iterator with location " + loc + " and rod name of " + ((RMDTrack)rod).getName());
            }
        }
    }
}


