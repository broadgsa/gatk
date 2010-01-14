package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
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
        SeekableRODIterator iterator = iteratorPool.iterator( new MappedStreamSegment(shard.getGenomeLoc()) );
        return iterator;
    }

    /**
     * Close the specified iterator, returning it to the pool.
     * @param iterator Iterator to close.
     */
    public void close( SeekableRODIterator iterator ) {
        this.iteratorPool.release(iterator);        
    }

}

/**
 * A pool of reference-ordered data iterators.
 */
class ReferenceOrderedDataPool extends ResourcePool<SeekableRODIterator,SeekableRODIterator> {
    private final ReferenceOrderedData<? extends ReferenceOrderedDatum> rod;

    public ReferenceOrderedDataPool( ReferenceOrderedData<? extends ReferenceOrderedDatum> rod ) {
        this.rod = rod;
    }

    /**
     * Create a new iterator from the existing reference-ordered data.  This new iterator is expected
     * to be completely independent of any other iterator.
     * @return The newly created resource.
     */
    public SeekableRODIterator createNewResource() {
        return rod.iterator();
    }

    /**
     * Finds the best existing ROD iterator from the pool.  In this case, the best existing ROD is defined as
     * the first one encountered that is at or before the given position.
     * @param segment @{inheritedDoc}
     * @param resources @{inheritedDoc}
     * @return @{inheritedDoc}
     */
    public SeekableRODIterator selectBestExistingResource( DataStreamSegment segment, List<SeekableRODIterator> resources ) {
        if( !(segment instanceof MappedStreamSegment) )
            throw new StingException("Reference-ordered data cannot utilitize unmapped segments.");

        GenomeLoc position = ((MappedStreamSegment)segment).locus;
        //#########################################
//## System.out.printf("Searching for iterator at locus %s; %d resources available%n", position, resources.size());
        for( SeekableRODIterator iterator: resources ) {
//##System.out.printf("Examining iterator at position %s [last query location: %s]%n", iterator.position(),iterator.lastQueryLocation());
            if( (iterator.position() == null && iterator.hasNext()) ||
                (iterator.position() != null && iterator.position().isBefore(position)) )
                return iterator;
        }
//##System.out.printf("Failed to find iterator at locus %s%n", position);
        return null;
    }

    /**
     * In this case, the iterator is the resource.  Pass it through.
     */
    public SeekableRODIterator createIteratorFromResource( DataStreamSegment segment, SeekableRODIterator resource ) {
        return resource;
    }

    /**
     * Don't worry about closing the resource; let the file handles expire naturally for the moment.
     */
    public void closeResource(  SeekableRODIterator resource ) {
        
    }
}


