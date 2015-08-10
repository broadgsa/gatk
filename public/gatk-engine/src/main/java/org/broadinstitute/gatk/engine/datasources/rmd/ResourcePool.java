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
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * A pool of open resources, all of which can create a closeable iterator.
 */
abstract class ResourcePool <T,I extends CloseableIterator> {
    /**
     * Sequence dictionary.
     */
    protected final SAMSequenceDictionary referenceSequenceDictionary;

    /**
     * Builder/parser for GenomeLocs.
     */
    protected final GenomeLocParser genomeLocParser;

    /**
     * All iterators of this reference-ordered data.
     */
    private List<T> allResources = new ArrayList<T>();

    /**
     * All iterators that are not currently in service.
     */
    private List<T> availableResources = new ArrayList<T>();

    /**
     * Which iterators are assigned to which pools.
     */
    private Map<I,T> resourceAssignments = new HashMap<I,T>();

    protected ResourcePool(SAMSequenceDictionary referenceSequenceDictionary,GenomeLocParser genomeLocParser) {
        this.referenceSequenceDictionary = referenceSequenceDictionary;
        this.genomeLocParser = genomeLocParser;
    }

    /**
     * Get an iterator whose position is before the specified location.  Create a new one if none exists.
     * @param segment Target position for the iterator.
     * @return An iterator that can traverse the selected region.  Should be able to iterate concurrently with other
     *         iterators from tihs pool.
     */
    public I iterator( DataStreamSegment segment ) {
        // Grab the first iterator in the list whose position is before the requested position.
        T selectedResource = null;
        synchronized (this) {
            selectedResource = selectBestExistingResource(segment, availableResources);

            // No iterator found?  Create another.  It is expected that
            // each iterator created will have its own file handle.
            if (selectedResource == null) {
                selectedResource = createNewResource();
                addNewResource(selectedResource);
            }

            // Remove the iterator from the list of available iterators.
            availableResources.remove(selectedResource);
        }


        I iterator = createIteratorFromResource(segment, selectedResource);

        // also protect the resource assignment
        synchronized (this) {
            // Make a note of this assignment for proper releasing later.
            resourceAssignments.put(iterator, selectedResource);
        }

        return iterator;
    }

    /**
     * Release the lock on the given iterator, returning it to the pool.
     * @param iterator Iterator to return to the pool.
     */
    public void release( I iterator ) {
        synchronized(this) {
            // Find and remove the resource from the list of allocated resources.
            T resource = resourceAssignments.get( iterator );
            Object obj = resourceAssignments.remove(iterator);

            // Close the iterator.
            iterator.close();

            // make sure we actually removed the assignment
            if (obj == null)
                    throw new ReviewedGATKException("Failed to remove resource assignment; target key had no associated value in the resource assignment map");
            // Return the resource to the pool.
            if( !allResources.contains(resource) )
                throw new ReviewedGATKException("Iterator does not belong to the given pool.");
            availableResources.add(resource);
        }
    }

    /**
     * Add a resource to the list of available resources.  Useful if derived classes
     * want to seed the pool with a set of at a given time (like at initialization).
     * @param resource The new resource to add.
     */
    protected void addNewResource( T resource ) {
        synchronized(this) {
            allResources.add(resource);
            availableResources.add(resource);
        }
    }

    /**
     * If no appropriate resources are found in the pool, the system can create a new resource.
     * Delegate the creation of the resource to the subclass.
     * @return The new resource created.
     */
    protected abstract T createNewResource();

    /**
     * Find the most appropriate resource to acquire the specified data.
     * @param segment The data over which the resource is required.
     * @param availableResources A list of candidate resources to evaluate.
     * @return The best choice of the availableResources, or null if no resource meets the criteria.
     */
    protected abstract T selectBestExistingResource( DataStreamSegment segment, List<T> availableResources );

    /**
     * Create an iterator over the specified resource.
     * @param position The bounds of iteration.  The first element of the iterator through the last element should all
     *                 be in the range described by position.
     * @param resource The resource from which to derive the iterator.
     * @return A new iterator over the given data.
     */
    protected abstract I createIteratorFromResource( DataStreamSegment position, T resource );

    /**
     * Retire this resource from service.
     * @param resource The resource to retire.
     */
    protected abstract void closeResource(T resource);

    /**
     * Operating stats...get the number of total iterators.  Package-protected
     * for unit testing.
     * @return An integer number of total iterators.
     */
    int numIterators() {
        return allResources.size();
    }

    /**
     * Operating stats...get the number of available iterators.  Package-protected
     * for unit testing.
     * @return An integer number of available iterators.
     */
    int numAvailableIterators() {
        return availableResources.size();
    }

}

