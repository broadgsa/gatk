package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.RODIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.ArrayList;
/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 10:55:26 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A pool of open iterators.  Currently highly specialized to RODs, but could theoretically be
 * generalized to a pool of arbitrary seekable, closeable iterators.  Not thread-safe.
 */
class IteratorPool {
    private final ReferenceOrderedData<? extends ReferenceOrderedDatum> rod;

    /**
     * All iterators of this reference-ordered data.
     */
    private List<RODIterator> allIterators = new ArrayList<RODIterator>();

    /**
     * All iterators that are not currently in service.
     */
    private List<RODIterator> availableIterators = new ArrayList<RODIterator>();

    /**
     * Create a new iterator pool given the current ROD.
     * @param rod Reference-ordered data.
     */
    public IteratorPool( ReferenceOrderedData<? extends ReferenceOrderedDatum> rod ) {
        this.rod = rod;
    }

    /**
     * Get an iterator whose position is before the specified location.  Create a new one if none exists.
     * @param position Target position for the iterator.
     * @return
     */
    public RODIterator iterator( GenomeLoc position ) {
        // Grab the first iterator in the list whose position is before the requested position.
        RODIterator selectedIterator = null;
        for( RODIterator iterator: availableIterators ) {
            if( (iterator.position() == null && iterator.hasNext()) ||
                (iterator.position() != null && iterator.position().isBefore(position)) ) {
                selectedIterator = iterator;
                break;
            }
        }

        // No iterator found?  Create another.  It is expected that
        // each iterator created will have its own file handle.
        if( selectedIterator == null ) {
            selectedIterator = rod.iterator();
            allIterators.add(selectedIterator);
        }

        // Remove the iterator from the list of available iterators.
        if( availableIterators.contains(selectedIterator) )
            availableIterators.remove(selectedIterator);

        return selectedIterator;
    }

    /**
     * Close the given iterator, returning it to the pool.
     * @param iterator Iterator to return to the pool.
     */
    public void close( RODIterator iterator ) {
        if( !allIterators.contains(iterator) )
            throw new StingException("Iterator does not belong to the given pool.");
        availableIterators.add(iterator);
    }

    /**
     * Operating stats...get the number of total iterators.  Package-protected
     * for unit testing.
     * @return An integer number of total iterators.
     */
    int numIterators() {
        return allIterators.size();
    }

    /**
     * Operating stats...get the number of available iterators.  Package-protected
     * for unit testing.
     * @return An integer number of available iterators.
     */
    int numAvailableIterators() {
        return availableIterators.size();
    }

}
