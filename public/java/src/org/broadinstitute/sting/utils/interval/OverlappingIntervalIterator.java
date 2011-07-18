/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Oct 7, 2010
 * Time: 2:40:02 PM
 * To change this template use File | Settings | File Templates.
 */

/** This class provides an adapter to Iterator<GenomeLoc> that returns only (parts of) underlying iterator's
 * intervals overlapping with specified "master set" of bounding intervals. The underlying iterator must return
 * NON-overlapping intervals in coordinate-sorted order, otherwise the behavior is unspecified. If the master set is represented by
 * another interval iterator, it should return sorted and NON-overlapping intervals.
 *
  */
public class OverlappingIntervalIterator implements Iterator<GenomeLoc> {
    PushbackIterator<GenomeLoc> iter = null;
    PushbackIterator<GenomeLoc> boundBy = null;

    GenomeLoc prefetchedOverlap = null;
    GenomeLoc currentBound = null;
    GenomeLoc currentInterval = null;

    
    /** Creates new overlapping iterator that will internally traverse <code>intervals</code> and return only
     * overlaps of those with set of intervals returned by <code>boundBy</code>.
      * @param intervals
     * @param boundBy
     */
    public OverlappingIntervalIterator(Iterator<GenomeLoc> intervals, Iterator<GenomeLoc> boundBy) {
        this.iter = new PushbackIterator<GenomeLoc>(intervals);
        this.boundBy = new PushbackIterator<GenomeLoc>(boundBy);

        if ( iter.hasNext() && boundBy.hasNext() ) {
            currentInterval = iter.next(); // load first interval
            currentBound = boundBy.next(); // load first bounding interval
            fetchNextOverlap();
        }
    }

    /** Traverses both iterators in sync, until the first overlap between the two is reached. If no overlap is found
     * until the end of the either of the two streams, leaves prefetchedOverlap set to null
     */
    private void fetchNextOverlap() {

        prefetchedOverlap = null;
 //       System.out.println("Fetching... (interval="+currentInterval+"; bound="+currentBound+")");
        while ( currentInterval != null && currentBound != null ) {

            if ( currentInterval.isBefore(currentBound) ) {
//                System.out.println(currentInterval +" is before "+currentBound );
                if ( ! iter.hasNext() ) currentInterval = null;
                else currentInterval = iter.next();
                continue;
            }

            if ( currentInterval.isPast(currentBound) ) {
//                System.out.println(currentInterval +" is past "+currentBound );
                if ( ! boundBy.hasNext() ) currentBound = null;
                else currentBound = boundBy.next();
                continue;
            }

            // we are at this point only if currentInterval overlaps with currentBound

            prefetchedOverlap = currentInterval.intersect(currentBound);
//            System.out.println("Fetched next overlap: "+prefetchedOverlap);
            // now we need to advance at least one of the iterators, so that we would not
            // call the same overlap again

            // however we still do not know if we are done with either current interval or current bound, because
            // two special situations are possible:
            //
            //     1) next interval overlaps with              2) current interval also overlaps with
            //        the same bounding interval;                  next bounding interval; note that
            //        note that in this case next                  in this case next bound necessarily
            //        interval necessarily starts before           starts before the next interval
            //        the next bound
            //
            //   curr. int     next int.                              curr. int
            //     -----        ------                        --------------------------
            //      -------------------                     ---------         -------------
            //           curr. bound                       curr. bound          next bound

            // To solve this issue we update either only currentInterval or only currentBound to their next value,
            // whichever of those next values (intervals) comes first on the reference genome;
            // the rest of the traversal to the next overlap will be performed on the next invocation of
            // fetchNextOverlap().

            advanceToNearest();

            break; // now that we computed the overlap and advanced (at least one of) the intervals/bounds to
                   // the next location, we are done - bail out from the loop.
        }

    }

    private void advanceToNearest() {
        if ( ! iter.hasNext() ) {
           currentBound = boundBy.hasNext() ? boundBy.next() : null;
        } else {
            if ( ! boundBy.hasNext() ) currentInterval = iter.hasNext() ? iter.next() : null;
            else {
                // both intervals and bounds have next value available; let's check which comes first:
                GenomeLoc nextInterval = iter.next();
                GenomeLoc nextBound = boundBy.next();

                if ( nextInterval.compareTo(nextBound) < 0 ) {
                    currentInterval = nextInterval;
                    boundBy.pushback(nextBound);
                } else {
                    currentBound = nextBound;
                    iter.pushback(nextInterval);
                }

            }
        }
    }

    /**
     * Returns <tt>true</tt> if the iteration has more elements. (In other
     * words, returns <tt>true</tt> if <tt>next</tt> would return an element
     * rather than throwing an exception.)
     *
     * @return <tt>true</tt> if the iterator has more elements.
     */
    public boolean hasNext() {
        return prefetchedOverlap != null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next element in the iteration.
     * @throws java.util.NoSuchElementException
     *          iteration has no more elements.
     */
    public GenomeLoc next() {
        if ( prefetchedOverlap == null )
            throw new java.util.NoSuchElementException("Illegal call to next(): Overlapping iterator has no more overlaps");
        GenomeLoc ret = prefetchedOverlap; // cache current prefetched overlap
        fetchNextOverlap(); // prefetch next overlap
        return ret ;
    }

    /**
     * Removes from the underlying collection the last element returned by the
     * iterator (optional operation).  This method can be called only once per
     * call to <tt>next</tt>.  The behavior of an iterator is unspecified if
     * the underlying collection is modified while the iteration is in
     * progress in any way other than by calling this method.
     *
     * @throws UnsupportedOperationException if the <tt>remove</tt>
     *                                       operation is not supported by this Iterator.
     * @throws IllegalStateException         if the <tt>next</tt> method has not
     *                                       yet been called, or the <tt>remove</tt> method has already
     *                                       been called after the last call to the <tt>next</tt>
     *                                       method.
     */
    public void remove() {
        throw new UnsupportedOperationException("remove() method is not supported by OverlappingIntervalIterator");
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
