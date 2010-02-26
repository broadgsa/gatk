package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.refdata.utils.FlashBackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.Iterator;
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
    private final ReferenceOrderedData rod;

    /**
     * A pool of iterators for navigating through the genome.
     */
    private ReferenceOrderedDataPool iteratorPool = null;

    /**
     * Create a new reference-ordered data source.
     * @param rod
     */
    public ReferenceOrderedDataSource( ReferenceOrderedData rod) {
        this.rod = rod;
        this.iteratorPool = new ReferenceOrderedDataPool( rod );
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
    public ReferenceOrderedData getReferenceOrderedData() {
        return this.rod;
    }

    /**
     * Seek to the specified position and return an iterator through the data.
     * @param shard Shard that points to the selected position.
     * @return Iterator through the data.
     */
    public Iterator seek( Shard shard ) {
        DataStreamSegment dataStreamSegment = shard.getGenomeLocs().size() != 0 ? new MappedStreamSegment(shard.getGenomeLocs().get(0)) : new EntireStream();
        FlashBackIterator RODIterator = iteratorPool.iterator(dataStreamSegment);
        return RODIterator;
    }

    /**
     * Close the specified iterator, returning it to the pool.
     * @param iterator Iterator to close.
     */
    public void close( FlashBackIterator iterator ) {
        this.iteratorPool.release(iterator);        
    }

}

/**
 * A pool of reference-ordered data iterators.
 */
class ReferenceOrderedDataPool extends ResourcePool<FlashBackIterator, FlashBackIterator> {
    private final ReferenceOrderedData<? extends ReferenceOrderedDatum> rod;
    public ReferenceOrderedDataPool( ReferenceOrderedData<? extends ReferenceOrderedDatum> rod ) {
        this.rod = rod;
    }

    /**
     * Create a new iterator from the existing reference-ordered data.  This new iterator is expected
     * to be completely independent of any other iterator.
     * @return The newly created resource.
     */
    public FlashBackIterator createNewResource() {
        return new FlashBackIterator(rod.iterator());
    }

    /**
     * Finds the best existing ROD iterator from the pool.  In this case, the best existing ROD is defined as
     * the first one encountered that is at or before the given position.
     * @param segment @{inheritedDoc}
     * @param resources @{inheritedDoc}
     * @return @{inheritedDoc}
     */
    public FlashBackIterator selectBestExistingResource( DataStreamSegment segment, List<FlashBackIterator> resources ) {
        if(segment instanceof MappedStreamSegment) {
            GenomeLoc position = ((MappedStreamSegment)segment).getFirstLocation();

            for( FlashBackIterator RODIterator : resources ) {

                if( (RODIterator.position() == null && RODIterator.hasNext()) ||
                    (RODIterator.position() != null && RODIterator.position().isBefore(position)) )
                    return RODIterator;
//                if (RODIterator.position() != null && RODIterator.canFlashBackTo(position)) {
//                    RODIterator.flashBackTo(position);
//                    return RODIterator;
//                }

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
    public FlashBackIterator createIteratorFromResource( DataStreamSegment segment, FlashBackIterator resource ) {
        return resource;
    }

    /**
     * kill the buffers in the iterator
     */
    public void closeResource( FlashBackIterator resource ) {
        resource.close();            
    }
}


