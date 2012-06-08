/*
 * Copyright (c) 2011 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.BitSetUtils;
import org.broadinstitute.sting.utils.clipping.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.BitSet;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements StandardCovariate {

    private int mismatchesContextSize;
    private int insertionsContextSize;
    private int deletionsContextSize;

    private byte LOW_QUAL_TAIL;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        insertionsContextSize = RAC.INSERTIONS_CONTEXT_SIZE;
        deletionsContextSize = RAC.DELETIONS_CONTEXT_SIZE;

        LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
        
        if (mismatchesContextSize <= 0 || insertionsContextSize <= 0 || deletionsContextSize <= 0)
            throw new UserException(String.format("Context Size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Insertions: %d Deletions:%d", mismatchesContextSize, insertionsContextSize, deletionsContextSize));

    }

    @Override
    public CovariateValues getValues(final GATKSAMRecord read) {
        int l = read.getReadLength();
        BitSet[] mismatches = new BitSet[l];
        BitSet[] insertions = new BitSet[l];
        BitSet[] deletions = new BitSet[l];

        GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);   // Write N's over the low quality tail of the reads to avoid adding them into the context
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        for (int i = 0; i < clippedRead.getReadLength(); i++) {
            mismatches[i] = contextWith(bases, i, mismatchesContextSize);
            insertions[i] = contextWith(bases, i, insertionsContextSize);
            deletions[i] = contextWith(bases, i, deletionsContextSize);
        }

        if (negativeStrand) {
            reverse(mismatches);
            reverse(insertions);
            reverse(deletions);
        }
        return new CovariateValues(mismatches, insertions, deletions);
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String keyFromBitSet(BitSet key) {
        if (key == null)    // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return BitSetUtils.dnaFrom(key);
    }

    @Override
    public BitSet bitSetFromKey(Object key) {
        return BitSetUtils.bitSetFrom((String) key);
    }

    @Override
    public int numberOfBits() {
        return Long.bitCount(-1L);
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param offset      the position in the read to calculate the context for
     * @param contextSize context size to use building the context
     * @return the bitSet representing the Context
     */
    private BitSet contextWith(byte[] bases, int offset, int contextSize) {
        BitSet result = null;
        if (offset - contextSize + 1 >= 0) {
            final byte[] context = Arrays.copyOfRange(bases, offset - contextSize + 1, offset + 1);
            if (!BaseUtils.containsBase(context, BaseUtils.N))
                result = BitSetUtils.bitSetFrom(context);
        }
        return result;
    }

    /**
     * Reverses the given array in place.
     *
     * @param array any array
     */
    private static void reverse(final Object[] array) {
        final int arrayLength = array.length;
        for (int l = 0, r = arrayLength - 1; l < r; l++, r--) {
            final Object temp = array[l];
            array[l] = array[r];
            array[r] = temp;
        }
    }
}
