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

package org.broadinstitute.gatk.utils;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Utilities for bitset conversion
 *
 * @author Mauricio Carneiro
 * @since 3/5/12
 */
public class BitSetUtils {

    static final private byte NBITS_LONG_REPRESENTATION = 64;                   // the number of bits used in the long version to represent the bit set (necessary for the two's complement representation of negative numbers)
    static final private byte NBITS_SHORT_REPRESENTATION = 16;                  // the number of bits used in the short version to represent the bit set (necessary for the two's complement representation of negative numbers)

    /**
     * Creates an long out of a bitset
     *
     * @param bitSet the bitset
     * @return a long from the bitset representation
     */
    public static long longFrom(final BitSet bitSet) {
        return longFrom(bitSet, NBITS_LONG_REPRESENTATION);
    }

    /**
     * Creates a short integer from a bitset
     *
     * @param bitSet the bitset
     * @return a short from the bitset representation
     */
    public static short shortFrom(final BitSet bitSet) {
        return (short) longFrom(bitSet, NBITS_SHORT_REPRESENTATION);
    }

    /**
     * Cretes an integer with any number of bits (up to 64 -- long precision) from a bitset
     *
     * @param bitSet the bitset
     * @param nBits  the number of bits to be used for this representation
     * @return an integer with nBits from the bitset representation
     */
    public static long longFrom(final BitSet bitSet, final int nBits) {
        long number = 0;
        for (int bitIndex = bitSet.nextSetBit(0); bitIndex >= 0 && bitIndex <= nBits; bitIndex = bitSet.nextSetBit(bitIndex + 1))
            number |= 1L << bitIndex;

        return number;
    }

    /**
     * Creates a BitSet representation of a given long
     *
     * @param number the number to turn into a bitset
     * @return a bitset representation of the long
     */
    public static BitSet bitSetFrom(long number) {
        return bitSetFrom(number, NBITS_LONG_REPRESENTATION);
    }

    /**
     * Creates a BitSet representation of a given short
     *
     * @param number the number to turn into a bitset
     * @return a bitset representation of the short
     */
    public static BitSet bitSetFrom(short number) {
        BitSet result = shortCache.get(number);
        if (result == null) {
            result = bitSetFrom(number, NBITS_SHORT_REPRESENTATION);
            shortCache.put(number, result);
        }
        return result;
    }
    // use a static cache for shorts (but not for longs, because there could be a lot of entries)
    private static final Map<Short, BitSet> shortCache = new HashMap<Short, BitSet>(2 * Short.MAX_VALUE);

    /**
     * Creates a BitSet representation of an arbitrary integer (number of bits capped at 64 -- long precision)
     *
     * @param number the number to turn into a bitset
     * @param nBits  the number of bits to use as precision for this conversion
     * @return a bitset representation of the integer
     */
    public static BitSet bitSetFrom(long number, int nBits) {
        BitSet bitSet = new BitSet(nBits);
        boolean isNegative = number < 0;
        int bitIndex = 0;
        while (number != 0) {
            if (number % 2 != 0)
                bitSet.set(bitIndex);
            bitIndex++;
            number /= 2;
        }
        if (isNegative) {
            boolean foundFirstSetBit = false;
            for (int i = bitSet.nextSetBit(0); i < nBits && i >= 0; i++) {
                boolean bit = bitSet.get(i);
                if (!foundFirstSetBit && bit)
                    foundFirstSetBit = true;    // maintain all bits until the first 1 is found (inclusive)
                else if (foundFirstSetBit)
                    bitSet.flip(i);             // flip every other bit up to NBITS_REPRESENTATION
            }
        }
        return bitSet;
    }
}
