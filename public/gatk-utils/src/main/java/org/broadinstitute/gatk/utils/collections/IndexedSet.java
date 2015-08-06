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

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.util.*;

/**
* Set set where each element can be reference by a unique integer index that runs from
*     0 to the size of the set - 1.
*
* @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
*/
public class IndexedSet<E> extends AbstractSet<E> implements Set<E> {

    /**
     * Elements stored in an array-list by their index.
     */
    private final ArrayList<E> elements;

    /**
     * A unmodifiable view to the element list. Initially {@code null} it is thread-unsafe lazy instantiated
     * when requested first time through {@link #asList}. Therefore typically it is shared by invoking code but
     * there could be some extra copies (rare though) in multi-thread runs.
     */
    private transient List<E> unmodifiableElementsListView;

    /**
     * Quick element to index lookup map.
     * <p>
     *  Uses a primitive int value map for efficiency sake.
     * </p>
     */
    private final Object2IntMap<E> indexByElement;

    /**
     * Creates an empty indexed set indicating the expected number of elements.
     *
     * @param initialCapacity the initial number of elements.
     */
    public IndexedSet(final int initialCapacity) {
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new Object2IntOpenHashMap<>(initialCapacity);
    }

    /**
     * Creates a new sample list from a existing collection of elements.
     *
     * <p>
     *     Elements will be indexed as they appear in the input array. Repeats will be ignored.
     * </p>
     *
     * @param values the original sample list.
     *
     * @throws IllegalArgumentException
     * if {@code values} array is {@code null} itself, or it contains {@code null}.
     */
    @SuppressWarnings("unchecked")
    public IndexedSet(final Collection<E> values) {
        if (values == null)
            throw new IllegalArgumentException("input values cannot be null");

        final int initialCapacity = values.size();
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new Object2IntOpenHashMap<>(initialCapacity);
        int nextIndex = 0;
        for (final E value : values) {
            if (value == null)
                throw new IllegalArgumentException("null element not allowed: index == " + nextIndex);
            if (indexByElement.containsKey(value))
                continue;
            indexByElement.put(value, nextIndex++);
            elements.add(value);
        }
    }

    /**
     * Creates a new sample list from a existing array of elements.
     *
     * <p>
     *     Elements will be indexed as they appear in the collection. Repeats will be ignored.
     * </p>
     *
     * @param values the original sample list.
     *
     * @throws IllegalArgumentException
     * if {@code values} collection is {@code null} itself, or it contains {@code null}.
     */
    @SuppressWarnings("unchecked")
    public IndexedSet(final E ... values) {
        if (values == null)
            throw new IllegalArgumentException("input values cannot be null");

        final int initialCapacity = values.length;
        elements = new ArrayList<>(initialCapacity);
        indexByElement = new Object2IntOpenHashMap<>(initialCapacity);
        int nextIndex = 0;
        for (final E value : values) {
            if (value == null)
                throw new IllegalArgumentException("null element not allowed: index == " + nextIndex);
            if (indexByElement.containsKey(value))
                continue;
            indexByElement.put(value, nextIndex++);
            elements.add(value);
        }
    }

    /**
     * Returns a list view of the elements in the set.
     *
     * <p>
     *     Elements are sorted by their index within the set.
     * </p>
     *
     * <p>
     *     This view changes as the indexed set changes but it cannot be used to update its contents.
     *     In such case a {@link UnsupportedOperationException} exception will be thrown if the calling
     *     code tries to tho just that.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<E> asList() {
        if (unmodifiableElementsListView == null)
            unmodifiableElementsListView = Collections.unmodifiableList(elements);
        return unmodifiableElementsListView;
    }

    /**
     * Throws an exception if an index is out of bounds.
     *
     * <p>
     *     An element index is valid iff is within [0,{@link #size()}).
     * </p>
     *
     * @param index the query index.
     *
     * @throws IllegalArgumentException {@code index} is out of bounds.
     */
    protected void checkIndex(final int index) {
        if (index < 0)
            throw new IllegalArgumentException("the index cannot be negative: " + index);
        if (index >= size())
            throw new IllegalArgumentException("the index is equal or larger than the list length: " + index + " >= " + size());
    }

    @Override
    public Iterator<E> iterator() {
        return asList().iterator();
    }

    /**
     * Returns number of elements in the set.
     * @return never {@code null}.
     */
    @Override
    public int size() {
        return elements.size();
    }

    /**
     *
     * @param o
     * @return {@code true} iff {@code o} is in
     */
    @Override
    @SuppressWarnings("all")
    public boolean contains(final Object o) {
        return o != null && indexByElement.containsKey(o);
    }

    /**
     * Adds a new element to the set.
     *
     * <p>
     *     If the element was already in th set nothing will happen and the method will return {@code false}. However,
     *     if the element is new to this set, it will assigned the next index available (equal to the size before addition).
     *     The method will return {@code true} in this case.
     * </p>
     *
     * @param o the object to add.
     *
     * @throw IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code true} iff the set was modified by this operation.
     */
    @Override
    public boolean add(final E o) {
        if (o == null)
            throw new IllegalArgumentException("the input argument cannot be null");
        if (contains(o))
            return false;
        final int nextIndex = size();
        elements.add(o);
        indexByElement.put(o, nextIndex);
        return true;
    }

    /**
     * Removes an element from the set.
     *
     * <p>
     *     If the element was not present in the set, nothing happens and the method return false. However,
     *     if the element is new to this set, it will be assigned the next index available (equal to the size
     *     before addition).
     *     The method will return {@code true} in this case.
     * </p>
     *
     * @param o the object to add.
     *
     * @throw IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code true} iff the set was modified by this operation.
     */   @Override
    public boolean remove(final Object o) {
        final int index = indexByElement.removeInt(o);
        if (index == -1)
            return false;
        elements.remove(index);
        indexByElement.remove(o);
        final ListIterator<E> it = elements.listIterator(index);
        int nextIndex = index;
        while (it.hasNext())
            indexByElement.put(it.next(),nextIndex++);
        return true;
    }

    /**
     * Removes all elements in the set.
     */
    @Override
    public void clear() {
        elements.clear();
        indexByElement.clear();
    }

    /**
     * Compares this with another indexed set.
     * @param o the other object to compare to.
     * @return {@code false} unless {@code o} is a indexed-set that contains the same elements in the same order.
     */
    @Override
    public boolean equals(final Object o) {
        if (o == this)
            return true;
        if (o == null)
            return false;
        if (!(o instanceof IndexedSet<?>))
            return false;

        final IndexedSet<?> other = (IndexedSet<?>)o;

        return equals(other);
    }

    /**
     * Compare to another indexed set.
     *
     * @param other the target indexed set.
     *
     * @throws java.lang.IllegalArgumentException if {@code other} is {@code null}.
     *
     * @return {@code true} iff {@other} is not {@code null}, and contains exactly the same elements
     * (as compared using {@link Object#equals} a this set with matching indices.
     */
    public boolean equals(final IndexedSet<?> other) {
        if (other == null)
            throw new IllegalArgumentException("other cannot be null");
        final ArrayList<?> otherElements = other.elements;

        final int elementCount = elements.size();
        if (otherElements.size() != elementCount)
            return false;
        for (int i = 0; i < elementCount; i++)
            if (!elements.get(i).equals(otherElements.get(i)))
                return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = 1;

        for (final E element : elements)
            result = 31 * result + (element == null ? 0 : element.hashCode());
        return result;
    }

    /**
     * Returns the element given its index within the set.
     * @param index the target element's index.
     *
     * @throws IllegalArgumentException if {@code index} is not valid; in [0,{@link #size()}).
     *
     * @return never {@code null}; as null is not a valid element.
     */
    public E get(final int index) {
        checkIndex(index);
        return elements.get(index);
    }

    /**
     * Returns the index of an object.
     * @param o the object of interest.
     *
     * @throws IllegalArgumentException if {@code o} is {@code null}.
     *
     * @return {@code -1} if such an object is not an element of this set, otherwise is index in the set thus a
     * values within [0,{@link #size()}).
     */
    public int indexOf(final E o) {
        if (o == null)
            throw new IllegalArgumentException("the query object cannot be null");
        return indexByElement.containsKey(o) ? indexByElement.getInt(o) : -1;
    }

}
