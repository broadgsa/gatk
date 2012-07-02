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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.collections;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: July 1, 2012
 */

public class IntegerIndexedNestedHashMap<T> {

    protected final Object[] data;

    protected final int numDimensions;
    protected final int[] indexFactors;

    public IntegerIndexedNestedHashMap(final int... dimensions) {
        numDimensions = dimensions.length;
        if ( numDimensions == 0 )
            throw new ReviewedStingException("There must be at least one dimension to an IntegerIndexedNestedHashMap");

        indexFactors = new int[numDimensions];
        indexFactors[0] = 1;
        int totalSize = dimensions[0];
        for ( int i = 1; i < numDimensions; i++ ) {
            indexFactors[i] = indexFactors[i-1] * dimensions[i-1];
            totalSize *= dimensions[i];
        }
        data = new Object[totalSize];
    }

    public T get(final int... keys) {
        return (T)data[calculateIndex(keys)];
    }

    public synchronized void put(final T value, final int... keys) { // WARNING! value comes before the keys!
        data[calculateIndex(keys)] = value;
    }

    protected int calculateIndex(final int... keys) {
        int index = keys[0];
        for ( int i = 1; i < numDimensions; i++ )
            index += keys[i] * indexFactors[i];
        return index;
    }

    public List<T> getAllValues() {
        final int dataLength = data.length;
        final List<T> result = new ArrayList<T>(dataLength);
        for ( int i = 0; i < dataLength; i++ ) {
            final Object value = data[i];
            if ( value != null )
                result.add((T)value);
        }
        return result;
    }

    public static class Leaf {
        public final int[] keys;
        public final Object value;

        public Leaf(final int[] keys, final Object value) {
            this.keys = keys;
            this.value = value;
        }
    }

    public List<Leaf> getAllLeaves() {
        final int dataLength = data.length;
        final List<Leaf> result = new ArrayList<Leaf>(dataLength);
        for ( int i = 0; i < dataLength; i++ ) {
            final Object value = data[i];
            if ( value != null )
                result.add(new Leaf(reverseEngineerIndex(i), value));
        }
        return result;
    }

    private int[] reverseEngineerIndex(int index) {
        final int[] keys = new int[numDimensions];
        for ( int i = numDimensions - 1; i >= 0; i-- ) {
            final int key = index / indexFactors[i];
            index -= (key * indexFactors[i]);
            keys[i] = key;
        }
        return keys;
    }
}
