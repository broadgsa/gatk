package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.iterators.PeekingIterator;
import org.broadinstitute.sting.gatk.iterators.PushbackIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;

import java.util.*;

public class MergingIterator<ROD extends ReferenceOrderedDatum> implements Iterator<RODRecordList<ROD>>, Iterable<RODRecordList<ROD>> {
    PriorityQueue<Element> queue = new PriorityQueue<Element>();

    private class Element implements Comparable<Element> {
        public SeekableRODIterator it = null;
        //public E value = null;
        public GenomeLoc nextLoc = null;

        public Element(Iterator<RODRecordList<ROD>> it) {
            if ( it instanceof SeekableRODIterator ) {
                this.it = (SeekableRODIterator)it;
                if ( ! it.hasNext() ) throw new StingException("Iterator is empty");
                update();
            } else {
                throw new StingException("Iterator passed to MergingIterator is not SeekableRODIterator");
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

        public RODRecordList<ROD> next() {
            RODRecordList<ROD> value = it.next();
            update();
            return value;
        }
    }

    public Iterator<RODRecordList<ROD>> iterator() {
        return this;
    }

    public MergingIterator() {
        ;
    }

    public MergingIterator(Iterator<RODRecordList<ROD>> it) {
         add(it);
    }

    public MergingIterator(Collection<Iterator<RODRecordList<ROD>>> its) {
        for ( Iterator<RODRecordList<ROD>> it : its ) {
            add(it);
        }
    }

    /** If the iterator is non-empty (hasNext() is true), put it into the queue. The next location the iterator
     * will be after a call to next() is peeked into and cached as queue's priority value.
     * @param it
     */
    public void add(Iterator<RODRecordList<ROD>> it) {
        if ( it.hasNext() )
            queue.add(new Element(it));
    }

    public boolean hasNext() {
        return ! queue.isEmpty();
    }

    public RODRecordList<ROD> next() {
        Element e = queue.poll();
        RODRecordList<ROD> value = e.next(); // next() will also update next location cached by the Element

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

    public Collection<RODRecordList<ROD>> allElementsLTE(RODRecordList<ROD> elt) {
        return allElementsLTE(elt, true);
    }

    public Collection<RODRecordList<ROD>> allElementsLTE(RODRecordList<ROD> elt, boolean includeElt) {
        LinkedList<RODRecordList<ROD>> all = new LinkedList<RODRecordList<ROD>>();

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
