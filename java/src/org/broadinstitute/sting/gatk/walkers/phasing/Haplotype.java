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

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Arrays;

public class Haplotype extends BaseArray implements Cloneable {
    public Haplotype(byte[] bases) {
        super(bases);
    }

    private Haplotype(Byte[] bases) {
        super(bases);
    }

    public Haplotype(Haplotype other) {
        super(other);
    }

    public Haplotype(BaseArray baseArr) {
        super(baseArr.bases);

        if (baseArr.getNonNullIndices().length != baseArr.bases.length)
            throw new ReviewedStingException("Should NEVER call Haplotype ctor with null bases!");
    }

    public void updateBase(int index, Byte base) {
        if (base == null) {
            throw new ReviewedStingException("Internal error: CANNOT have null for a missing Haplotype base!");
        }
        super.updateBase(index, base);
    }

    public Haplotype clone() {
        try {
            super.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return new Haplotype(this);
    }

    // Returns a new Haplotype containing the portion of this Haplotype between the specified fromIndex, inclusive, and toIndex, exclusive.

    public Haplotype subHaplotype(int fromIndex, int toIndex) {
        return new Haplotype(Arrays.copyOfRange(bases, fromIndex, Math.min(toIndex, size())));
    }
}
