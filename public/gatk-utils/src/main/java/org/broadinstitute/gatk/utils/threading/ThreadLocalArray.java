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

package org.broadinstitute.gatk.utils.threading;

import java.lang.reflect.Array;

/**
 * ThreadLocal implementation for arrays
 *
 * Example usage:
 *
 * private ThreadLocal<byte[]> threadLocalByteArray = new ThreadLocalArray<byte[]>(length, byte.class);
 * ....
 * byte[] byteArray = threadLocalByteArray.get();
 *
 * @param <T> the type of the array itself (eg., int[], double[], etc.)
 *
 * @author David Roazen
 */
public class ThreadLocalArray<T> extends ThreadLocal<T> {
    private int arraySize;
    private Class arrayElementType;

    /**
     * Create a new ThreadLocalArray
     *
     * @param arraySize desired length of the array
     * @param arrayElementType type of the elements within the array (eg., Byte.class, Integer.class, etc.)
     */
    public ThreadLocalArray( int arraySize, Class arrayElementType ) {
        super();

        this.arraySize = arraySize;
        this.arrayElementType = arrayElementType;
    }

    @Override
    @SuppressWarnings("unchecked")
    protected T initialValue() {
        return (T)Array.newInstance(arrayElementType, arraySize);
    }
}
