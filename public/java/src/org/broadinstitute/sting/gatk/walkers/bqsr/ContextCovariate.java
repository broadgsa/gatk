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
import org.broadinstitute.sting.utils.clipping.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipping.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements StandardCovariate {

    private int mismatchesContextSize;
    private int indelsContextSize;

    // the maximum context size (number of bases) permitted; we need to keep the leftmost base free so that values are
    // not negative and we reserve 4 more bits to represent the length of the context; it takes 2 bits to encode one base.
    static final private int MAX_DNA_CONTEXT = 13;
    private byte LOW_QUAL_TAIL;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        indelsContextSize = RAC.INDELS_CONTEXT_SIZE;
        if (mismatchesContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        if (indelsContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("indels_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, indelsContextSize));

        LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
        
        if (mismatchesContextSize <= 0 || indelsContextSize <= 0)
            throw new UserException(String.format("Context size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Indels: %d", mismatchesContextSize, indelsContextSize));
    }

    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {

        // TODO -- wrong: fix me
        final GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);   // Write N's over the low quality tail of the reads to avoid adding them into the context
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        final int readLength = clippedRead.getReadLength();
        for (int i = 0; i < readLength; i++) {
            final int indelKey = contextWith(bases, i, indelsContextSize);
            values.addCovariate(contextWith(bases, i, mismatchesContextSize), indelKey, indelKey, (negativeStrand ? readLength - i - 1 : i));
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String formatKey(final int key) {
        if (key == -1)    // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return contextFromKey(key);
    }

    @Override
    public int keyFromValue(final Object value) {
        return keyFromContext((String) value);
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param offset      the position in the read to calculate the context for
     * @param contextSize context size to use building the context
     * @return the key representing the context
     */
    private int contextWith(final byte[] bases, final int offset, final int contextSize) {
        final int start = offset - contextSize + 1;
        return (start >= 0) ? keyFromContext(bases, start, offset + 1) : -1;
    }

    public static int keyFromContext(final String dna) {
        return keyFromContext(dna.getBytes(), 0, dna.length());
    }

    /**
     * Creates a int representation of a given dna string.
     *
     * @param dna    the dna sequence
     * @param start  the start position in the byte array (inclusive)
     * @param end    the end position in the array (exclusive)
     * @return the key representing the dna sequence
     */
    public static int keyFromContext(final byte[] dna, final int start, final int end) {

        // TODO -- bit fiddle to ge this all working in a single call to the method (mask out length, shift, OR length back in)

        int key = end - start;
        int bitOffset = 4;
        for (int i = start; i < end; i++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
            if (baseIndex == -1)                    // ignore non-ACGT bases
                return -1;
            key |= (baseIndex << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    /**
     * Converts a key into the dna string representation.
     *
     * @param key    the key representing the dna sequence
     * @return the dna sequence represented by the key
     */
    public static String contextFromKey(final int key) {
        if (key < 0)
            throw new ReviewedStingException("dna conversion cannot handle negative numbers. Possible overflow?");

        final int length = key & 15;               // the first 4 bits represent the length (in bp) of the context
        int mask = 48;                             // use the mask to pull out bases
        int offset = 4;

        StringBuilder dna = new StringBuilder();
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask = mask << 2;                      // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }
}
