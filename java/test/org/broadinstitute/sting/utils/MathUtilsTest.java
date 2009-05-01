package org.broadinstitute.sting.utils;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;


import java.io.File;

/**
 * Basic unit test for MathUtils
 */
public class MathUtilsTest extends BaseTest {
    @BeforeClass
    public static void init() { }

    /**
     * Tests that we get the right values from the binomial distribution
     */
    @Test
    public void testBinomialProbability() {
        logger.warn("Executing testBinomialProbability");

        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(2, 3, 0.5), 0.375, 0.0001) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(10, 100, 0.5), 1.365543e-17, 1e-18) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(73, 217, 0.02), 4.521904e-67, 1e-68) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(100, 300, 0.02), 9.27097e-91, 1e-92) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(150, 300, 0.98), 6.462892e-168, 1e-169) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(120, 300, 0.98), 3.090054e-221, 1e-222) == 0);
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.binomialProbability(112, 300, 0.98), 2.34763e-236, 1e-237) == 0);
    }

    /**
     * Tests that we get the right values from the multinomial distribution
     */
    @Test
    public void testMultinomialProbability() {
        logger.warn("Executing testMultinomialProbability");

        int[] counts0 = { 2, 0, 1 };
        double[] probs0 = { 0.33, 0.33, 0.34 };
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.multinomialProbability(counts0, probs0), 0.111078, 1e-6) == 0);

        int[] counts1 = { 10, 20, 30 };
        double[] probs1 = { 0.25, 0.25, 0.50 };
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.multinomialProbability(counts1, probs1), 0.002870301, 1e-9) == 0);

        int[] counts2 = { 38, 82, 50, 36 };
        double[] probs2 = { 0.25, 0.25, 0.25, 0.25 };
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.multinomialProbability(counts2, probs2), 1.88221e-09, 1e-10) == 0);

        int[] counts3 = { 1, 600, 1 };
        double[] probs3 = { 0.33, 0.33, 0.34 };
        Assert.assertTrue(MathUtils.compareDoubles(MathUtils.multinomialProbability(counts3, probs3), 5.20988e-285, 1e-286) == 0);
    }
}