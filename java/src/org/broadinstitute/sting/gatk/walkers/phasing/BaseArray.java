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

import org.broadinstitute.sting.utils.BaseUtils;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public abstract class BaseArray implements Comparable<BaseArray> {
    protected Byte[] bases;

    public BaseArray(byte[] bases) {
        this.bases = new Byte[bases.length];
        for (int i = 0; i < bases.length; i++)
            this.bases[i] = bases[i];
    }

    public BaseArray(Byte[] bases) {
        this.bases = Arrays.copyOf(bases, bases.length);
    }

    public BaseArray(int length) {
        this.bases = new Byte[length];
        Arrays.fill(this.bases, null);
    }

    public BaseArray(BaseArray other) {
        this(other.bases);
    }

    public void updateBase(int index, Byte base) {
        bases[index] = base;
    }

    public Byte getBase(int index) {
        return bases[index];
    }

    public int size() {
        return bases.length;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(bases.length);
        for (Byte b : bases)
            sb.append(b != null ? (char) b.byteValue() : "_");

        return sb.toString();
    }

    public int compareTo(BaseArray that) {
        int sz = this.bases.length;
        if (sz != that.bases.length)
            return (sz - that.bases.length);

        for (int i = 0; i < sz; i++) {
            Byte thisBase = this.getBase(i);
            Byte thatBase = that.getBase(i);
            if (thisBase == null || thatBase == null) {
                if (thisBase == null && thatBase != null) {
                    return -1;
                }
                else if (thisBase != null && thatBase == null) {
                    return 1;
                }
            }
            else if (!BaseUtils.extendedBasesAreEqual(thisBase, thatBase)) {
                return thisBase - thatBase;
            }
        }
        return 0;
    }

    public int[] getNonNullIndices() {
        List<Integer> nonNull = new LinkedList<Integer>();
        for (int i = 0; i < bases.length; i++) {
            if (getBase(i) != null)
                nonNull.add(i);
        }

        int[] nonNullArray = new int[nonNull.size()];
        int index = 0;
        for (int i : nonNull)
            nonNullArray[index++] = i;
        return nonNullArray;
    }
}
