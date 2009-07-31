package org.broadinstitute.sting.utils;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MathUtils;


import java.io.File;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Basic unit test for MathUtils
 */
public class ListUtilsTest extends BaseTest {
    @BeforeClass
    public static void init() { }

    /**
     * Tests that the random index selection is working correctly
     */
    @Test
    public void testRandomIndicesWithReplacement() {
        logger.warn("Executing testRandomIndicesWithReplacement");

        // Check that the size of the list returned is correct
        Assert.assertTrue(ListUtils.sampleIndicesWithReplacement(0, 5).size() == 0);
        Assert.assertTrue(ListUtils.sampleIndicesWithReplacement(1, 5).size() == 1);
        Assert.assertTrue(ListUtils.sampleIndicesWithReplacement(5, 5).size() == 5);
        Assert.assertTrue(ListUtils.sampleIndicesWithReplacement(1000, 5).size() == 1000);

        // Check that the list contains only the k element range that as asked for - no more, no less
        List<Integer> Five = new ArrayList<Integer>();
        Collections.addAll(Five, 0, 1, 2, 3, 4);
        List<Integer> BigFive = ListUtils.sampleIndicesWithReplacement(10000, 5);
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
        List<Integer> BigFive = ListUtils.sampleIndicesWithReplacement(10000, 5);
        List<Character> BigFiveAlpha = ListUtils.sliceListByIndices(BigFive, FiveAlpha);
        Assert.assertTrue(BigFiveAlpha.containsAll(FiveAlpha));
        Assert.assertTrue(FiveAlpha.containsAll(BigFiveAlpha));
    }
}