/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.phasing;

import java.util.NoSuchElementException;

/*
  DoublyLinkedList class is a doubly-linked list, which allows O(1) traversal to next and previous elements in the list.
  It is UNIQUE in the fact that its iterator (BidirectionalIterator) can be cloned
  to save the current pointer for a later time (while the original iterator can continue to iterate).
 */
public class CloneableIteratorLinkedList<E> {
    private CloneableIteratorDoublyLinkedNode<E> first;
    private CloneableIteratorDoublyLinkedNode<E> last;
    private int size;

    public CloneableIteratorLinkedList() {
        this.first = null;
        this.last = null;
        this.size = 0;
    }

    public boolean isEmpty() {
        return first == null;
    }

    public int size() {
        return size;
    }

    public void addFirst(E e) {
        CloneableIteratorDoublyLinkedNode<E> newNode = new CloneableIteratorDoublyLinkedNode<E>(e);

        if (isEmpty())
            last = newNode;
        else {
            first.previous = newNode;
            newNode.next = first;
        }
        first = newNode;

        size++;
    }

    public void addLast(E e) {
        CloneableIteratorDoublyLinkedNode<E> newNode = new CloneableIteratorDoublyLinkedNode<E>(e);

        if (isEmpty())
            first = newNode;
        else {
            last.next = newNode;
            newNode.previous = last;
        }
        last = newNode;

        size++;
    }

    public E removeFirst() {
        if (isEmpty())
            throw new NoSuchElementException();
        E e = first.element;

        if (first.next == null)
            last = null;
        else
            first.next.previous = null;
        first = first.next;

        size--;
        return e;
    }

    public E removeLast() {
        if (isEmpty())
            throw new NoSuchElementException();
        E e = last.element;

        if (last.previous == null)
            first = null;
        else
            last.previous.next = null;
        last = last.previous;

        size--;
        return e;
    }

    public E getFirst() {
        if (isEmpty())
            throw new NoSuchElementException();

        return first.element;
    }

    public E getLast() {
        if (isEmpty())
            throw new NoSuchElementException();

        return last.element;
    }

    public E peek() {
        if (isEmpty())
            return null;

        return getFirst();
    }

    public E remove() {
        return removeFirst();
    }

    public boolean add(E e) {
    	addLast(e);
        return true;
    }

    public CloneableIterator<E> iterator() {
        return new CloneableIterator<E>(this);
    }


    private static class CloneableIteratorDoublyLinkedNode<E> {
        private E element = null;
        private CloneableIteratorDoublyLinkedNode<E> next = null;
        private CloneableIteratorDoublyLinkedNode<E> previous = null;

        public CloneableIteratorDoublyLinkedNode(E element) {
            this.element = element;
            this.next = null;
            this.previous = null;
        }
    }


    /*
    This iterator is unique since it can be cloned to save the current pointer for a later time (while the original iterator can continue to iterate).
     */
    public static class CloneableIterator<E> implements Cloneable {
        private CloneableIteratorDoublyLinkedNode<E> nextNode;
        private CloneableIteratorDoublyLinkedNode<E> lastNode;

        private CloneableIterator(CloneableIteratorDoublyLinkedNode<E> nextNode, CloneableIteratorDoublyLinkedNode<E> lastNode) {
            this.nextNode = nextNode;
            this.lastNode = lastNode;
        }

        private CloneableIterator(CloneableIteratorLinkedList<E> list) {
            this(list.first, list.last);
        }

        public boolean hasNext() {
            return nextNode != null;
        }

        public E next() {
            if (!hasNext())
                throw new NoSuchElementException();

            E e = nextNode.element;
            nextNode = nextNode.next;
            return e;
        }

        public boolean hasPrevious() {
            if (nextNode != null)
                return nextNode.previous != null;

            return lastNode != null;
        }

        public E previous() {
            if (!hasPrevious())
                throw new NoSuchElementException();

            if (nextNode != null)
                nextNode = nextNode.previous;
            else
                nextNode = lastNode;

            return nextNode.element;
        }

        public CloneableIterator<E> clone() {
            try {
                super.clone();
            } catch (CloneNotSupportedException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            return new CloneableIterator<E>(nextNode, lastNode);
        }
    }
}
