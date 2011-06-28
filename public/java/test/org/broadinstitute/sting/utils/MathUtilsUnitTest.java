/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils;


import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import org.broadinstitute.sting.BaseTest;


import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Basic unit test for MathUtils
 */
public class MathUtilsUnitTest extends BaseTest {
    @BeforeClass
    public void init() { }

    /**
     * Tests that we get the right values from the binomial distribution
     */
    @Test
    public void testBinomialProbability() {
        logger.warn("Executing testBinomialProbability");

        Assert.assertEquals(MathUtils.binomialProbability(3, 2, 0.5), 0.375, 0.0001);
        Assert.assertEquals(MathUtils.binomialProbability(100, 10, 0.5), 1.365543e-17, 1e-18);
        Assert.assertEquals(MathUtils.binomialProbability(217, 73, 0.02), 4.521904e-67, 1e-68);
        Assert.assertEquals(MathUtils.binomialProbability(300, 100, 0.02), 9.27097e-91, 1e-92);
        Assert.assertEquals(MathUtils.binomialProbability(300, 150, 0.98), 6.462892e-168, 1e-169);
        Assert.assertEquals(MathUtils.binomialProbability(300, 120, 0.98), 3.090054e-221, 1e-222);
        Assert.assertEquals(MathUtils.binomialProbability(300, 112, 0.98), 2.34763e-236, 1e-237);
    }

    /**
     * Tests that we get the right values from the multinomial distribution
     */
    @Test
    public void testMultinomialProbability() {
        logger.warn("Executing testMultinomialProbability");

        int[] counts0 = { 2, 0, 1 };
        double[] probs0 = { 0.33, 0.33, 0.34 };
        Assert.assertEquals(MathUtils.multinomialProbability(counts0, probs0), 0.111078, 1e-6);

        int[] counts1 = { 10, 20, 30 };
        double[] probs1 = { 0.25, 0.25, 0.50 };
        Assert.assertEquals(MathUtils.multinomialProbability(counts1, probs1), 0.002870301, 1e-9);

        int[] counts2 = { 38, 82, 50, 36 };
        double[] probs2 = { 0.25, 0.25, 0.25, 0.25 };
        Assert.assertEquals(MathUtils.multinomialProbability(counts2, probs2), 1.88221e-09, 1e-10);

        int[] counts3 = { 1, 600, 1 };
        double[] probs3 = { 0.33, 0.33, 0.34 };
        Assert.assertEquals(MathUtils.multinomialProbability(counts3, probs3), 5.20988e-285, 1e-286);
    }

    /**
     * Tests that the random index selection is working correctly
     */
    @Test
    public void testRandomIndicesWithReplacement() {
        logger.warn("Executing testRandomIndicesWithReplacement");

        // Check that the size of the list returned is correct
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 0).size() == 0);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 1).size() == 1);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 5).size() == 5);
        Assert.assertTrue(MathUtils.sampleIndicesWithReplacement(5, 1000).size() == 1000);

        // Check that the list contains only the k element range that as asked for - no more, no less
        List<Integer> Five = new ArrayList<Integer>();
        Collections.addAll(Five, 0, 1, 2, 3, 4);
        List<Integer> BigFive = MathUtils.sampleIndicesWithReplacement(5, 10000);
        Assert.assertTrue(BigFive.containsAll(Five));
        Assert.assertTrue(Five.containsAll(BigFive));
    }

    /**
     * Tests that we get the right values from the multinomial distribution
     */
    @Test
    public void testSliceListByIndices() {
        logger.warn("Executing testSliceListByIndices");

        // Check that the list contains only the k element range that as asked for - no more, no less but now
        // use the index list to pull elements from another list using sliceListByIndices
        List<Integer> Five = new ArrayList<Integer>();
        Collections.addAll(Five, 0, 1, 2, 3, 4);
        List<Character> FiveAlpha = new ArrayList<Character>();
        Collections.addAll(FiveAlpha, 'a', 'b', 'c', 'd', 'e');
        List<Integer> BigFive = MathUtils.sampleIndicesWithReplacement(5, 10000);
        List<Character> BigFiveAlpha = MathUtils.sliceListByIndices(BigFive, FiveAlpha);
        Assert.assertTrue(BigFiveAlpha.containsAll(FiveAlpha));
        Assert.assertTrue(FiveAlpha.containsAll(BigFiveAlpha));
    }

    /** Tests that we correctly compute mean and standard deviation from a stream of numbers */
    @Test
    public void testRunningAverage() {
        logger.warn("Executing testRunningAverage");

        int [] numbers = {1,2,4,5,3,128,25678,-24};
        MathUtils.RunningAverage r = new MathUtils.RunningAverage();

        for ( int i = 0 ; i < numbers.length ; i++ ) r.add((double)numbers[i]);

        Assert.assertEquals((long)numbers.length, r.observationCount());
        Assert.assertTrue(r.mean()- 3224.625 < 2e-10 );
        Assert.assertTrue(r.stddev()-9072.6515881128 < 2e-10);
    }

    @Test
    public void testLog10Gamma() {
        logger.warn("Executing testLog10Gamma");

        Assert.assertEquals(MathUtils.log10Gamma(4.0), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10), 5.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Gamma(10654), 38280.53, 1e-2);
    }

    @Test
    public void testLog10BinomialCoefficient() {
        logger.warn("Executing testLog10BinomialCoefficient");

        Assert.assertEquals(MathUtils.log10BinomialCoefficient(4, 2), 0.7781513, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(10, 3), 2.079181, 1e-6);
        Assert.assertEquals(MathUtils.log10BinomialCoefficient(103928, 119), 400.2156, 1e-4);
    }

    @Test
    public void testFactorial() {
        logger.warn("Executing testFactorial");
        Assert.assertEquals((int) MathUtils.factorial(4), 24);
        Assert.assertEquals((int) MathUtils.factorial(10), 3628800);
        Assert.assertEquals((int) MathUtils.factorial(12), 479001600);
    }

    @Test
    public void testLog10Factorial() {
        logger.warn("Executing testLog10Factorial");
        Assert.assertEquals(MathUtils.log10Factorial(4), 1.380211, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(10), 6.559763, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(12), 8.680337, 1e-6);
        Assert.assertEquals(MathUtils.log10Factorial(200), 374.8969, 1e-3);
        Assert.assertEquals(MathUtils.log10Factorial(12342), 45138.26, 1e-1);
    }

}
