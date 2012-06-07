/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.utils.codecs.bcf2;

import com.google.java.contract.Requires;

import java.util.EnumSet;

/**
 * BCF2 types and associated information
 *
 * @author depristo
 * @since 05/12
 */
public enum BCF2Type {
    INT8 (1, 1, 0xFFFFFF80,        -127,        127), // todo -- confirm range
    INT16(2, 2, 0xFFFF8000,      -32767,      32767),
    INT32(3, 4, 0x80000000, -2147483647, 2147483647),
    FLOAT(5, 4, 0x7F800001),
    CHAR (7);

    private final int id;
    private final Object missingJavaValue;
    private final int missingBytes;
    private final int sizeInBytes;
    private final long minValue, maxValue;

    BCF2Type(final int id) {
        this(id, -1, 0, 0, 0);
    }

    BCF2Type(final int id, final int sizeInBytes, final int missingBytes) {
        this(id, sizeInBytes, missingBytes, 0, 0);
    }

    BCF2Type(final int id, final int sizeInBytes, final int missingBytes, final long minValue, final long maxValue) {
        this.id = id;
        this.sizeInBytes = sizeInBytes;
        this.missingJavaValue = null;
        this.missingBytes = missingBytes;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    /**
     * How many bytes are used to represent this type on disk?
     * @return
     */
    public int getSizeInBytes() {
        return sizeInBytes;
    }

    /**
     * The ID according to the BCF2 specification
     * @return
     */
    public int getID() { return id; }

    /**
     * Can we encode value v in this type, according to its declared range.
     *
     * Only makes sense for integer values
     *
     * @param v
     * @return
     */
    @Requires("INTEGERS.contains(this)")
    public final boolean withinRange(final long v) { return v >= minValue && v <= maxValue; }

    /**
     * Return the java object (aka null) that is used to represent a missing value for this
     * type in Java
     *
     * @return
     */
    public Object getMissingJavaValue() { return missingJavaValue; }

    /**
     * The bytes (encoded as an int) that are used to represent a missing value
     * for this type in BCF2
     *
     * @return
     */
    public int getMissingBytes() { return missingBytes; }

    /**
     * An enum set of the types that might represent Integer values
     */
    public final static EnumSet<BCF2Type> INTEGERS = EnumSet.of(INT8, INT16, INT32);
}
