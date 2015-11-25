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

package org.broadinstitute.gatk.utils.nanoScheduler;

import org.broadinstitute.gatk.utils.collections.ExpandingArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 12/19/12
 * Time: 3:53 PM
 *
 * This class makes some critical assumptions.  First is that the jobID of the first
 * job is 0.  If this isn't true the MapResultsQueue will certainly fail.
 */
public class MapResultsQueue<MapType> {
    //private final static boolean DEBUG = false;
    //private final static Logger logger = Logger.getLogger(MapResultsQueue.class);

    /**
     * Although naturally stored as priority blocking queue, this is actually quite expensive
     * due to the O(n log n) sorting calculation.  Since we know that the job ids start
     * at 0 and increment by 1 in each successive job, we store an array instead.  The
     * array is indexed by jobID, and contains the MapResult for that job id.  Because elements
     * can be added to the queue in any order, we need to use an expanding array list to
     * store the elements.
     */
    final ExpandingArrayList<MapResult<MapType>> queue = new ExpandingArrayList<MapResult<MapType>>(10000);

    /**
     * The jobID of the last job we've seen
     */
    int prevJobID = -1; // no jobs observed

    /**
     * Put mapResult into this MapResultsQueue, associated with its jobID
     * @param mapResult a non-null map result
     */
    public synchronized void put(final MapResult<MapType> mapResult) {
        if ( mapResult == null ) throw new IllegalArgumentException("mapResult cannot be null");

        // make sure that nothing is at the job id for map
        assert queue.size() < mapResult.getJobID() || queue.get(mapResult.getJobID()) == null;

        queue.set(mapResult.getJobID(), mapResult);
    }

    /**
     * Should we reduce the next value in the mapResultQueue?
     *
     * @return true if we should reduce
     */
    public synchronized boolean nextValueIsAvailable() {
        final MapResult<MapType> nextMapResult = queue.get(nextJobID());

        if ( nextMapResult == null ) {
            // natural case -- the next job hasn't had a value added yet
            return false;
        } else if ( nextMapResult.getJobID() != nextJobID() ) {
            // sanity check -- the job id at next isn't the one we expect
            throw new IllegalStateException("Next job ID " + nextMapResult.getJobID() + " is not == previous job id " + prevJobID + " + 1");
        } else {
            // there's a value at the next job id, so return true
            return true;
        }
    }

    /**
     * Get the next job ID'd be expect to see given our previous job id
     * @return the next job id we'd fetch to reduce
     */
    private int nextJobID() {
        return prevJobID + 1;
    }

    /**
     * Can only be called when nextValueIsAvailable is true
     * @return
     * @throws InterruptedException
     */
    // TODO -- does this have to be synchronized? -- I think the answer is no
    public synchronized MapResult<MapType> take() throws InterruptedException {
        final MapResult<MapType> result = queue.get(nextJobID());

        // make sure the value we've fetched has the right id
        assert result.getJobID() == nextJobID();

        prevJobID = result.getJobID();
        queue.set(prevJobID, null);

        return result;
    }
}
