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

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.Random;

/**
 * @author Mauricio Carneiro
 * @since 3/5/12
 */

public class BitSetUtilsUnitTest {
    private static int RANDOM_NUMBERS_TO_TRY = 87380;
    private static Random random;

    @BeforeClass
    public void init() {
        random = Utils.getRandomGenerator();
    }

    @Test(enabled = true)
    public void testLongBitSet() {
        long[] numbers = {0L, 1L, 428L, 65536L, 239847L, 4611686018427387903L, Long.MAX_VALUE, Long.MIN_VALUE, -1L, -2L, -7L, -128L, -65536L, -100000L};
        for (long n : numbers)
            Assert.assertEquals(BitSetUtils.longFrom(BitSetUtils.bitSetFrom(n)), n);

        for (int i = 0; i < RANDOM_NUMBERS_TO_TRY; i++) {
            long n = random.nextLong();
            Assert.assertEquals(BitSetUtils.longFrom(BitSetUtils.bitSetFrom(n)), n);    // Because class Random uses a seed with only 48 bits, this algorithm will not return all possible long values.
        }
    }

    @Test(enabled = true)
    public void testShortBitSet() {
        short[] numbers = {0, 1, 428, 25934, 23847, 16168, Short.MAX_VALUE, Short.MIN_VALUE, -1, -2, -7, -128, -12312, -31432};
        for (long n : numbers)
            Assert.assertEquals(BitSetUtils.shortFrom(BitSetUtils.bitSetFrom(n)), n);

        for (int i = 0; i < RANDOM_NUMBERS_TO_TRY; i++) {
            short n = (short) random.nextInt();
            Assert.assertEquals(BitSetUtils.shortFrom(BitSetUtils.bitSetFrom(n)), n);
        }
    }

    @Test(enabled = false)
    public void testDNAAndBitSetConversion() {
        String[] dna = {"AGGTGTTGT", "CCCCCCCCCCCCCC", "GGGGGGGGGGGGGG", "TTTTTTTTTTTTTT", "GTAGACCGATCTCAGCTAGT", "AACGTCAATGCAGTCAAGTCAGACGTGGGTT", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"};

        // Test all contexts of size 1-8.
        //for (long n = 0; n < RANDOM_NUMBERS_TO_TRY; n++)
        //    Assert.assertEquals(BitSetUtils.longFrom(BitSetUtils.bitSetFrom(ContextCovariate.contextFromKey(BitSetUtils.bitSetFrom(n)))), n);

        // Test the special cases listed in the dna array
        //for (String d : dna)
        //    Assert.assertEquals(BitSetUtils.dnaFrom(BitSetUtils.bitSetFrom(d)), d);
    }
}
