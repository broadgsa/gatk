package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RODIterator;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
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
    private final ReferenceOrderedData<? extends ReferenceOrderedDatum> rod;

    /**
     * A pool of iterators for navigating through the genome.
     */
    private IteratorPool iteratorPool = null;

    /**
     * Create a new reference-ordered data source.
     * @param rod
     */
    public ReferenceOrderedDataSource( ReferenceOrderedData<? extends ReferenceOrderedDatum> rod) {
        this.rod = rod;
        this.iteratorPool = new IteratorPool( rod );
    }

    /**
     * Return the name of the underlying reference-ordered data.
     * @return Name of the underlying rod.
     */
    public String getName() {
        return this.rod.getName();
    }

    /**
     * Seek to the specified position and return an iterator through the data.
     * @param shard Shard that points to the selected position.
     * @return Iterator through the data.
     */
    public Iterator seek( Shard shard ) {
        RODIterator iterator = iteratorPool.iterator(shard.getGenomeLoc());
        return iterator;
    }

    /**
     * Close the specified iterator, returning it to the pool.
     * @param iterator Iterator to close.
     */
    public void close( RODIterator iterator ) {
        this.iteratorPool.close(iterator);        
    }
}


