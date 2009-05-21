package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.BaseTest;
import org.junit.Test;
import org.junit.BeforeClass;
import org.junit.Assert;

public class BaseUtilsTest extends BaseTest {
    @BeforeClass
    public static void init() { }

    @Test
    public void testMostFrequentBaseFraction() {
        logger.warn("Executing testMostFrequentBaseFraction");
        
        compareFrequentBaseFractionToExpected("AAAAA", 1.0);
        compareFrequentBaseFractionToExpected("ACCG", 0.5);
        compareFrequentBaseFractionToExpected("ACCCCTTTTG", 4.0/10.0);
    }

    private void compareFrequentBaseFractionToExpected(String sequence, double expected) {
        double fraction = BaseUtils.mostFrequentBaseFraction(sequence.getBytes());
        Assert.assertTrue(MathUtils.compareDoubles(fraction, expected) == 0);
    }
}
