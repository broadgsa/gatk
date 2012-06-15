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
    private int insertionsContextSize;
    private int deletionsContextSize;

    private byte LOW_QUAL_TAIL;

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        mismatchesContextSize = RAC.MISMATCHES_CONTEXT_SIZE;
        insertionsContextSize = RAC.INSERTIONS_CONTEXT_SIZE;
        deletionsContextSize = RAC.DELETIONS_CONTEXT_SIZE;
        if (mismatchesContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("mismatches_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, mismatchesContextSize));
        if (insertionsContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("insertions_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, insertionsContextSize));
        if (deletionsContextSize > MAX_DNA_CONTEXT)
            throw new UserException.BadArgumentValue("deletions_context_size", String.format("context size cannot be bigger than %d, but was %d", MAX_DNA_CONTEXT, deletionsContextSize));

        LOW_QUAL_TAIL = RAC.LOW_QUAL_TAIL;
        
        if (mismatchesContextSize <= 0 || insertionsContextSize <= 0 || deletionsContextSize <= 0)
            throw new UserException(String.format("Context Size must be positive, if you don't want to use the context covariate, just turn it off instead. Mismatches: %d Insertions: %d Deletions:%d", mismatchesContextSize, insertionsContextSize, deletionsContextSize));
    }

    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {

        final GATKSAMRecord clippedRead = ReadClipper.clipLowQualEnds(read, LOW_QUAL_TAIL, ClippingRepresentation.WRITE_NS);   // Write N's over the low quality tail of the reads to avoid adding them into the context
        
        final boolean negativeStrand = clippedRead.getReadNegativeStrandFlag();
        byte[] bases = clippedRead.getReadBases();
        if (negativeStrand)
            bases = BaseUtils.simpleReverseComplement(bases);

        final int readLength = clippedRead.getReadLength();
        for (int i = 0; i < readLength; i++) {
            values.addCovariate(contextWith(bases, i, mismatchesContextSize), contextWith(bases, i, insertionsContextSize), contextWith(bases, i, deletionsContextSize), (negativeStrand ? readLength - i - 1 : i));
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return str;
    }

    @Override
    public String formatKey(final long key) {
        if (key == -1)    // this can only happen in test routines because we do not propagate null keys to the csv file
            return null;

        return contextFromKey(key);
    }

    @Override
    public long longFromKey(Object key) {
        return keyFromContext((String) key);
    }

    @Override
    public int numberOfBits() {
        return Integer.bitCount(Integer.MAX_VALUE);
    }

    /**
     * calculates the context of a base independent of the covariate mode (mismatch, insertion or deletion)
     *
     * @param bases       the bases in the read to build the context from
     * @param offset      the position in the read to calculate the context for
     * @param contextSize context size to use building the context
     * @return the key representing the context
     */
    private long contextWith(final byte[] bases, final int offset, final int contextSize) {
        final int start = offset - contextSize + 1;
        final long result;
        if (start >= 0)
            result = keyFromContext(bases, start, offset + 1);
        else
            result = -1L;
        return result;
    }

    public static long keyFromContext(final String dna) {
        return keyFromContext(dna.getBytes(), 0, dna.length());
    }

    /**
     * Creates a long representation of a given dna string.
     *
     * Warning: This conversion is limited to long precision, therefore the dna sequence cannot
     * be longer than 31 bases.
     *
     * The bit representation of a dna string is the simple:
     * 0 A      4 AA     8 CA
     * 1 C      5 AC     ...
     * 2 G      6 AG     1343 TTGGT
     * 3 T      7 AT     1364 TTTTT
     *
     * To convert from dna to number, we convert the dna string to base10 and add all combinations that
     * preceded the string (with smaller lengths).
     *
     * @param dna the dna sequence
     * @return the key representing the dna sequence
     */
    public static long keyFromContext(final byte[] dna, final int start, final int end) {
        final long preContext = combinationsPerLength[end - start - 1];      // the sum of all combinations that preceded the length of the dna string
        long baseTen = 0L;                                                   // the number in base_10 that we are going to use to generate the bit set
        for (int i = start; i < end; i++) {
            baseTen = (baseTen << 2);               // multiply by 4
            final int baseIndex = BaseUtils.simpleBaseToBaseIndex(dna[i]);
            if (baseIndex == -1)                    // ignore non-ACGT bases
                return -1L;
            baseTen += (long)baseIndex;
        }
        return baseTen + preContext;                // the number representing this DNA string is the base_10 representation plus all combinations that preceded this string length.
    }

    static final private int MAX_DNA_CONTEXT = 31;                              // the maximum context size (number of bases) permitted in the "long bitset" implementation of the DNA <=> BitSet conversion.
    static final long[] combinationsPerLength = new long[MAX_DNA_CONTEXT + 1];  // keeps the memoized table with the number of combinations for each given DNA context length
    static {
        for (int i = 0; i < MAX_DNA_CONTEXT + 1; i++)
            computeCombinationsFor(i);
    }

    /**
     * The sum of all combinations of a context of a given length from length = 0 to length.
     *
     * Memoized implementation of sum(4^i) , where i=[0,length]
     *
     * @param length the length of the DNA context
     */
    private static void computeCombinationsFor(final int length) {
        long combinations = 0L;
        for (int i = 1; i <= length; i++)
            combinations += (1L << 2 * i);        // add all combinations with 4^i ( 4^i is the same as 2^(2*i) )
        combinationsPerLength[length] = combinations;
    }

    /**
     * Converts a key into the dna string representation.
     *
     * Warning: This conversion is limited to long precision, therefore the dna sequence cannot
     * be longer than 31 bases.
     *
     * We calculate the length of the resulting DNA sequence by looking at the sum(4^i) that exceeds the
     * base_10 representation of the sequence. This is important for us to know how to bring the number
     * to a quasi-canonical base_4 representation, and to fill in leading A's (since A's are represented
     * as 0's and leading 0's are omitted).
     *
     * quasi-canonical because A is represented by a 0, therefore,
     * instead of : 0, 1, 2, 3, 10, 11, 12, ...
     * we have    : 0, 1, 2, 3, 00, 01, 02, ...
     *
     * but we can correctly decode it because we know the final length.
     *
     * @param key    the key representing the dna sequence
     * @return the dna sequence represented by the key
     */
    public static String contextFromKey(long key) {
        if (key < 0)
            throw new ReviewedStingException("dna conversion cannot handle negative numbers. Possible overflow?");

        final int length = contextLengthFor(key);  // the length of the context (the number of combinations is memoized, so costs zero to separate this into two method calls)
        key -= combinationsPerLength[length - 1];  // subtract the the number of combinations of the preceding context from the number to get to the quasi-canonical representation

        StringBuilder dna = new StringBuilder();
        while (key > 0) {                         // perform a simple base_10 to base_4 conversion (quasi-canonical)
            final byte base = (byte) (key & 3);   // equivalent to (key % 4)
            dna.append((char)BaseUtils.baseIndexToSimpleBase(base));
            key = key >> 2;     // divide by 4
        }
        for (int j = dna.length(); j < length; j++)
            dna.append('A');                          // add leading A's as necessary (due to the "quasi" canonical status, see description above)

        return dna.reverse().toString();              // make sure to reverse the string since we should have been pre-pending all along
    }

    /**
     * Calculates the length of the DNA context for a given base 10 number
     *
     * It is important to know the length given the base 10 number to calculate the number of combinations
     * and to disambiguate the "quasi-canonical" state.
     *
     * This method also calculates the number of combinations as a by-product, but since it memoizes the
     * results, a subsequent call to combinationsFor(length) is O(1).
     *
     * @param number the base 10 representation of the key
     * @return the length of the DNA context represented by this number
     */
    private static int contextLengthFor(final long number) {
        int length = 1;                                     // the calculated length of the DNA sequence given the base_10 representation of its BitSet.
        long combinations = combinationsPerLength[length];  // the next context (we advance it so we know which one was preceding it).
        while (combinations <= number) {                    // find the length of the dna string (length)
            length++;
            combinations = combinationsPerLength[length];   // calculate the next context
        }
        return length;
    }
}
