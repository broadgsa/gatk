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

/**
 * BCF2 types and information
 *
 * @author depristo
 * @since 05/12
 */
public enum BCFType {
    RESERVED_0,
    INT8(1, BCF2Constants.INT8_MISSING_VALUE, -127, 127), // todo -- confirm range
    INT16(2, BCF2Constants.INT16_MISSING_VALUE, -32767, 32767),
    INT32(4, BCF2Constants.INT32_MISSING_VALUE, -2147483647, 2147483647),
    RESERVED_4,
    FLOAT(4, BCF2Constants.FLOAT_MISSING_VALUE),
    RESERVED_6,
    CHAR;

    private final Object missingJavaValue;
    private final int missingBytes;
    private final int sizeInBytes;
    private final long minValue, maxValue;

    BCFType() {
        this(-1, 0, 0, 0);
    }

    BCFType(final int sizeInBytes, final int missingBytes) {
        this(sizeInBytes, missingBytes, 0, 0);
    }

    BCFType(final int sizeInBytes, final int missingBytes, final long minValue, final long maxValue) {
        this.sizeInBytes = sizeInBytes;
        this.missingJavaValue = null;
        this.missingBytes = missingBytes;
        this.minValue = minValue;
        this.maxValue = maxValue;
    }

    public int getSizeInBytes() {
        return sizeInBytes;
    }
    public int getID() { return ordinal(); }
    public final boolean withinRange(final long v) { return v >= minValue && v <= maxValue; }
    public Object getMissingJavaValue() { return missingJavaValue; }
    public int getMissingBytes() { return missingBytes; }
}
