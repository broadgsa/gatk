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

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 9/26/11
 */

public class ContextCovariate implements StandardCovariate {

    private int mismatchesContextSize;
    private int indelsContextSize;

    private int mismatchesKeyMask;
    private int indelsKeyMask;

    private static final int LENGTH_BITS = 4;
    private static final int LENGTH_MASK = 15;

    // temporary lists to use for creating context covariate keys
    private final ArrayList<Integer> mismatchKeys = new ArrayList<Integer>(200);
    private final ArrayList<Integer> indelKeys = new ArrayList<Integer>(200);

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

        mismatchesKeyMask = createMask(mismatchesContextSize);
        indelsKeyMask = createMask(indelsContextSize);
    }

    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {

        // store the original bases and then write Ns over low quality ones
        final byte[] originalBases = read.getReadBases();
        final GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);   // Write N's over the low quality tail of the reads to avoid adding them into the context
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        mismatchKeys.clear();
        indelKeys.clear();
        contextWith(bases, mismatchesContextSize, mismatchKeys, mismatchesKeyMask);
        contextWith(bases, indelsContextSize, indelKeys, indelsKeyMask);

        final int readLength = bases.length;
        for (int i = 0; i < readLength; i++) {
            final int indelKey = indelKeys.get(i);
            values.addCovariate(mismatchKeys.get(i), indelKey, indelKey, (negativeStrand ? readLength - i - 1 : i));
        }

        // put the original bases back in
        read.setReadBases(originalBases);
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

    private static int createMask(final int contextSize) {
        int mask = 0;
        // create 2*contextSize worth of bits
        for (int i = 0; i < contextSize; i++)
            mask = (mask << 2) | 3;
        // shift 4 bits to mask out the bits used to encode the length
        return mask << LENGTH_BITS;
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param contextSize context size to use building the context
     * @param keys        list to store the keys
     * @param mask        mask for pulling out just the context bits
     */
    private static void contextWith(final byte[] bases, final int contextSize, final ArrayList<Integer> keys, final int mask) {

        // the first contextSize-1 bases will not have enough previous context
        for (int i = 1; i < contextSize && i <= bases.length; i++)
            keys.add(-1);

        if (bases.length < contextSize)
            return;

        final int newBaseOffset = 2 * (contextSize - 1) + LENGTH_BITS;

        // get (and add) the key for the context starting at the first base
        int currentKey = keyFromContext(bases, 0, contextSize);
        keys.add(currentKey);

        // if the first key was -1 then there was an N in the context; figure out how many more consecutive contexts it affects
        int currentNPenalty = 0;
        if (currentKey == -1) {
            currentKey = 0;
            currentNPenalty = contextSize - 1;
            int offset = newBaseOffset;
            while (bases[currentNPenalty] != 'N') {
                final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentNPenalty]);
                currentKey |= (baseIndex << offset);
                offset -= 2;
                currentNPenalty--;
            }
        }

        final int readLength = bases.length;
        for (int currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(bases[currentIndex]);
            if (baseIndex == -1) {                    // ignore non-ACGT bases
                currentNPenalty = contextSize;
                currentKey = 0;                       // reset the key
            } else {
                // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
                currentKey = (currentKey >> 2) & mask;
                currentKey |= (baseIndex << newBaseOffset);
                currentKey |= contextSize;
            }

            if (currentNPenalty == 0) {
                keys.add(currentKey);
            } else {
                currentNPenalty--;
                keys.add(-1);
            }
        }
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
    private static int keyFromContext(final byte[] dna, final int start, final int end) {

        int key = end - start;
        int bitOffset = LENGTH_BITS;
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

        final int length = key & LENGTH_MASK;               // the first bits represent the length (in bp) of the context
        int mask = 48;                                      // use the mask to pull out bases
        int offset = LENGTH_BITS;

        StringBuilder dna = new StringBuilder();
        for (int i = 0; i < length; i++) {
            final int baseIndex = (key & mask) >> offset;
            dna.append((char)BaseUtils.baseIndexToSimpleBase(baseIndex));
            mask = mask << 2;                      // move the mask over to the next 2 bits
            offset += 2;
        }

        return dna.toString();
    }

    @Override
    public int maximumKeyValue() {
        // the maximum value is T (11 in binary) for each base in the context
        int length = Math.max(mismatchesContextSize, indelsContextSize);  // the length of the context
        int key = length;
        int bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }
}
