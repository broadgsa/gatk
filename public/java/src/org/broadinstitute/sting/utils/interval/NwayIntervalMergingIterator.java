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

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.PriorityQueue;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Oct 28, 2010
 * Time: 12:06:23 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * An adapter over a collection of underlying Iterator<GenomeLoc> objects (a single underlying iterator is allowed). Each
 * individual underlying iterator must serve its intervals in coordinate-sorted order or an exception will be thrown.
 * Intervals from individual underlying streams (iterators) are 1) merged into a single ordered stream; 2) each group of
 * overlapping intervals from that merged stream are merged into a single interval; each call to next() returns such
 * merged interval guaranteed to have no overlaps with the previous or next interval. 
 *
 */
public class NwayIntervalMergingIterator implements Iterator<GenomeLoc>, Iterable<GenomeLoc> {

    private PriorityQueue<Element> queue = null;
    private IntervalMergingRule myRule;

    public NwayIntervalMergingIterator(IntervalMergingRule rule) {
        myRule = rule;
        queue = new PriorityQueue<Element>();
    }

    public void add(Iterator<GenomeLoc> it) {
        Element e = new Element(it);
        if ( ! e.isEmpty() ) queue.add(e);
    }

    public Iterator<GenomeLoc> iterator() {
        return this;
    }

    /**
     * Returns <tt>true</tt> if the iteration has more elements. (In other
     * words, returns <tt>true</tt> if <tt>next</tt> would return an element
     * rather than throwing an exception.)
     *
     * @return <tt>true</tt> if the iterator has more elements.
     */
    public boolean hasNext() {
        return ! queue.isEmpty();  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next element in the iteration.
     * @throws java.util.NoSuchElementException
     *          iteration has no more elements.
     */
    public GenomeLoc next() {
        Element e = queue.poll();
        GenomeLoc result = e.current;

        // advance element (i.e. its underlying iterator) and reinsert into the queue
        e.advance();
        if ( ! e.isEmpty() ) queue.add(e);

        while ( ! queue.isEmpty () ) {
            e = queue.peek();

            if (result.overlapsP(e.current) ||  myRule == IntervalMergingRule.ALL && result.contiguousP(e.current)) {
                // we need to merge:
                result = result.merge(e.current);

                // remove current head of the queue that we just merged into the result:
                e = queue.poll();
                // advance element we just merged into the result and reinsert it into the queue (if it has any data left):
                e.advance();
                if ( ! e.isEmpty() ) queue.add(e);

            } else {
                // next element does not overlap with current result; we are done: return the result and that
                // next element will wait for next call to next()
                break;
            }

        }
        return result;  //To change body of implemented methods use File | Settings | File Templates.
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
        throw new UnsupportedOperationException("remove() method not supported by this iterator");
    }

    private class Element implements Comparable<Element> {
        private Iterator<GenomeLoc> it;
        private GenomeLoc current = null;

        private void advance() {
            if ( it.hasNext() ) {
                GenomeLoc next = it.next();
                if ( next.isBefore(current) ) {
                    throw new UserException("Interval list provided by underlying iterator "+it.getClass().getName() +" is out of order");
                }
                current = next;
            }
            else current = null;
        }

        public boolean isEmpty() { return current == null; }

        public Element(Iterator<GenomeLoc> it) {
            this.it = it;
            if ( this.it.hasNext() ) current = this.it.next();
        }

        /**
         * Compares this object with the specified object for order.  Returns a
         * negative integer, zero, or a positive integer as this object is less
         * than, equal to, or greater than the specified object.
         * <p/>
         * <p>The implementor must ensure <tt>sgn(x.compareTo(y)) ==
         * -sgn(y.compareTo(x))</tt> for all <tt>x</tt> and <tt>y</tt>.  (This
         * implies that <tt>x.compareTo(y)</tt> must throw an exception iff
         * <tt>y.compareTo(x)</tt> throws an exception.)
         * <p/>
         * <p>The implementor must also ensure that the relation is transitive:
         * <tt>(x.compareTo(y)&gt;0 &amp;&amp; y.compareTo(z)&gt;0)</tt> implies
         * <tt>x.compareTo(z)&gt;0</tt>.
         * <p/>
         * <p>Finally, the implementor must ensure that <tt>x.compareTo(y)==0</tt>
         * implies that <tt>sgn(x.compareTo(z)) == sgn(y.compareTo(z))</tt>, for
         * all <tt>z</tt>.
         * <p/>
         * <p>It is strongly recommended, but <i>not</i> strictly required that
         * <tt>(x.compareTo(y)==0) == (x.equals(y))</tt>.  Generally speaking, any
         * class that implements the <tt>Comparable</tt> interface and violates
         * this condition should clearly indicate this fact.  The recommended
         * language is "Note: this class has a natural ordering that is
         * inconsistent with equals."
         * <p/>
         * <p>In the foregoing description, the notation
         * <tt>sgn(</tt><i>expression</i><tt>)</tt> designates the mathematical
         * <i>signum</i> function, which is defined to return one of <tt>-1</tt>,
         * <tt>0</tt>, or <tt>1</tt> according to whether the value of
         * <i>expression</i> is negative, zero or positive.
         *
         * @param o the object to be compared.
         * @return a negative integer, zero, or a positive integer as this object
         *         is less than, equal to, or greater than the specified object.
         * @throws ClassCastException if the specified object's type prevents it
         *                            from being compared to this object.
         */
        public int compareTo(Element o) {
            if ( current == null ) return 1;
            if ( o.current == null ) return -1;
            return current.compareTo(o.current);  //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
