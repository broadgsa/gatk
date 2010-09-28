package org.broadinstitute.sting.utils;

import java.util.NoSuchElementException;

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

public class DoublyLinkedList<E> {
    private DoublyLinkedNode<E> first;
    private DoublyLinkedNode<E> last;
    private int size;

    public DoublyLinkedList() {
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
        DoublyLinkedNode<E> newNode = new DoublyLinkedNode<E>(e);

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
        DoublyLinkedNode<E> newNode = new DoublyLinkedNode<E>(e);

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

    public BidirectionalIterator<E> iterator() {
        return new BidirectionalIterator<E>(this);
    }


    private static class DoublyLinkedNode<E> {
        private E element = null;
        private DoublyLinkedNode<E> next = null;
        private DoublyLinkedNode<E> previous = null;

        public DoublyLinkedNode(E element) {
            this.element = element;
            this.next = null;
            this.previous = null;
        }
    }


    public static class BidirectionalIterator<E> implements Cloneable {
        private DoublyLinkedNode<E> nextNode;
        private DoublyLinkedNode<E> lastNode;

        private BidirectionalIterator(DoublyLinkedNode<E> nextNode, DoublyLinkedNode<E> lastNode) {
            this.nextNode = nextNode;
            this.lastNode = lastNode;
        }

        private BidirectionalIterator(DoublyLinkedList<E> list) {
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

        public BidirectionalIterator<E> clone() {
            try {
                super.clone();
            } catch (CloneNotSupportedException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            return new BidirectionalIterator<E>(nextNode, lastNode);
        }
    }
}
