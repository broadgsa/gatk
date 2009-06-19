package org.broadinstitute.sting.utils;

import java.util.ArrayList;
import java.util.Collection;

public class ExpandingArrayList<E> extends ArrayList<E> {
    public ExpandingArrayList() { super(); }
    public ExpandingArrayList(Collection<? extends E> c) { super(c); }
    public ExpandingArrayList(int initialCapacity) { super(initialCapacity); }

    /**
     * Returns the element at the specified position in this list.  If index > size,
     * returns null.  Otherwise tries to access the array
     * @param index
     * @return
     * @throws IndexOutOfBoundsException in index < 0
     */
    public E get(int index) throws IndexOutOfBoundsException {
        if ( index < size() )
            return super.get(index);
        else
            return null;
    }

    public E set(int index, E element) {
        if ( index >= size() ) {
            // We need to add null items until we can safely set index to element
            for ( int i = size(); i <= index; i++ )
                add(null);
        }

        return super.set(index, element);
    }
}