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

package org.broadinstitute.gatk.utils.collections;

import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.PriorityQueue;

public class RODMergingIterator implements Iterator<RODRecordList>, Iterable<RODRecordList> {
    PriorityQueue<Element> queue = new PriorityQueue<Element>();

    private class Element implements Comparable<Element> {
        public LocationAwareSeekableRODIterator it = null;
        public GenomeLoc nextLoc = null;

        public Element(Iterator<RODRecordList> it) {
            if ( it instanceof LocationAwareSeekableRODIterator) {
                this.it = (LocationAwareSeekableRODIterator)it;
                if ( ! it.hasNext() ) throw new ReviewedGATKException("Iterator is empty");
                update();
            } else {
                throw new ReviewedGATKException("Iterator passed to RODMergingIterator is not LocationAwareSeekableRODIterator");
            }
        }

        public Element update() {
 //           E prev = value;
            nextLoc = it.peekNextLocation(); // will return null if there is no next location
            return this;
        }

        public int compareTo(Element other) {
            if ( nextLoc == null ) {
                if ( other.nextLoc != null ) return 1; // null means no more data available, so its after any non-null position
                return 0;
            }
            if ( other.nextLoc == null ) return -1; // we can get to this point only if this.nextLoc != null

            return nextLoc.compareTo(other.nextLoc);
        }

        public RODRecordList next() {
            RODRecordList value = it.next();
            update();
            return value;
        }
    }

    public Iterator<RODRecordList> iterator() {
        return this;
    }

    public RODMergingIterator() {
        ;
    }

    public RODMergingIterator(Iterator<RODRecordList> it) {
         add(it);
    }

    public RODMergingIterator(Collection<Iterator<RODRecordList>> its) {
        for ( Iterator<RODRecordList> it : its ) {
            add(it);
        }
    }

    /** If the iterator is non-empty (hasNext() is true), put it into the queue. The next location the iterator
     * will be after a call to next() is peeked into and cached as queue's priority value.
     * @param it
     */
    public void add(Iterator<RODRecordList> it) {
        if ( it.hasNext() )
            queue.add(new Element(it));
    }

    public boolean hasNext() {
        return ! queue.isEmpty();
    }

    public RODRecordList next() {
        Element e = queue.poll();
        RODRecordList value = e.next(); // next() will also update next location cached by the Element

        if ( e.nextLoc != null ) // we have more data in the track
            queue.add(e); // add the element back to queue (note: its next location, on which priority is based, was updated

        //System.out.printf("Element is %s%n", e.value);
        return value;
    }

    /** Peeks into the genomic location of the record this iterator will return next.
     *
     * @return
     */
    public GenomeLoc peekLocation() {
        return queue.peek().nextLoc;
    }

    public Collection<RODRecordList> allElementsLTE(RODRecordList elt) {
        return allElementsLTE(elt, true);
    }

    public Collection<RODRecordList> allElementsLTE(RODRecordList elt, boolean includeElt) {
        LinkedList<RODRecordList> all = new LinkedList<RODRecordList>();

        if ( includeElt ) all.add(elt);
        
        while ( hasNext() ) {
            Element x = queue.peek();
            //System.out.printf("elt.compareTo(x) == %d%n", elt.compareTo(x));
            //System.out.printf("In allElementLTE%n");
            int cmp = elt.getLocation().compareTo(x.nextLoc);
            //System.out.printf("x=%s%n  elt=%s%n  => elt.compareTo(x) == %d%n", x, elt, cmp);
            if ( cmp >= 0 ) {
                //System.out.printf("  Adding element x=%s, size = %d%n", x, all.size());
                all.add(next());
                //System.out.printf("  Added size = %d%n", all.size());
            }
            else {
                //System.out.printf("breaking...%n");
                break;
            }
        }

        return all;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
