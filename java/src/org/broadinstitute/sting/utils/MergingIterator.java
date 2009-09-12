package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.iterators.PeekingIterator;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;

import java.util.*;

public class MergingIterator<E extends Comparable<E>> implements Iterator<E>, PeekingIterator<E>, Iterable<E> {
    PriorityQueue<Element> queue = new PriorityQueue<Element>();

    private class Element implements Comparable<Element> {
        public Iterator<E> it = null;
        public E value = null;

        public Element(Iterator<E> it) {
            this.it = it;
            update();
        }

        public Element update() {
            if ( ! it.hasNext() )
                throw new RuntimeException("it is empty");

            E prev = value;
            value = it.next();
            //System.out.printf("Updating %s to prev=%s, next=%s%n", this, prev, value);
            return this;
        }

        public int compareTo(Element other) {
            return value.compareTo(other.value);
        }
    }

    public Iterator<E> iterator() {
        return this;
    }

    public MergingIterator() {
        ;
    }

    public MergingIterator(Iterator<E> it) {
         add(it);
    }

    public MergingIterator(Collection<Iterator<E>> its) {
        for ( Iterator<E> it : its ) {
            add(it);
        }
    }

    public void add(Iterator<E> it) {
        if ( it.hasNext() )
            queue.add(new Element(it));
    }

    public boolean hasNext() {
        return ! queue.isEmpty();
    }

    public E next() {
        Element e = queue.poll();
        E value = e.value;

        if ( e.it != null && e.it.hasNext() )
            queue.add(new Element(e.it));

        //System.out.printf("Element is %s%n", e.value);
        return value;
    }

    public E peek() {
        return queue.peek().value;
    }

    public Collection<E> allElementsLTE(E elt) {
        return allElementsLTE(elt, true);
    }

    public Collection<E> allElementsLTE(E elt, boolean includeElt) {
        LinkedList<E> all = new LinkedList<E>();

        if ( includeElt ) all.add(elt);
        
        while ( hasNext() ) {
            E x = peek();
            //System.out.printf("elt.compareTo(x) == %d%n", elt.compareTo(x));
            //System.out.printf("In allElementLTE%n");
            int cmp = elt.compareTo(x);
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
