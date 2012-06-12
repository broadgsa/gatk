package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.walkers.bqsr.BQSRKeyManager;
import org.broadinstitute.sting.gatk.walkers.bqsr.ContextCovariate;
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
        random = GenomeAnalysisEngine.getRandomGenerator();
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

    @Test(enabled = true)
    public void testNumberOfBitsToRepresent() {
        Assert.assertEquals(BQSRKeyManager.numberOfBitsToRepresent(0), 0); // Make sure 0 elements need 0 bits to be represented
        Assert.assertEquals(BQSRKeyManager.numberOfBitsToRepresent(1), 1); // Make sure 1 element needs 1 bit to be represented
        Assert.assertEquals(BQSRKeyManager.numberOfBitsToRepresent(3), 2); // Make sure 3 elements need 2 bit to be represented

        for (int i = 1; i < 63; i++) {                                  // Can't test i == 63 because n1 is a negative number
            long n1 = 1L << i;
            long n2 = Math.abs(random.nextLong()) % n1;
            long n3 = n1 | n2;
            Assert.assertEquals(BQSRKeyManager.numberOfBitsToRepresent(n3), (n3 == n1) ? i : i + 1);
            Assert.assertEquals(BQSRKeyManager.numberOfBitsToRepresent(n1), i);
        }
    }
}
