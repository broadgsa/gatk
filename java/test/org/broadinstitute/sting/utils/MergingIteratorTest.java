// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.

import org.junit.Assert;
import org.junit.Test;
import org.junit.Before;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

/**
 * Basic unit test for RecalData
 */
public class MergingIteratorTest extends BaseTest {
    List<Integer> one_to_ten = Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    List<Integer> evens = Arrays.asList(2, 4, 6, 8, 10);
    List<Integer> three_four_five = Arrays.asList(3, 4, 5);
    List<Integer> bigs = Arrays.asList(100, 200, 500);

    MergingIterator<Integer> x = null;

    @Before
    public void before() {
        x = null;
    }

    @Test
    public void testOneToOne() {
        logger.warn("Executing testOneToOne");

        x = new MergingIterator<Integer>(one_to_ten.iterator());

        for ( int i = 1; i < 11; i++ ) {
            int elt = x.next();
            logger.warn(String.format("Testing %d == %d", i, elt));
            Assert.assertEquals(elt, i);
        }
        Assert.assertFalse(x.hasNext());
    }

    @Test
    public void testMerging10Evens() {
        logger.warn("Executing testMerging10Evens");

        x = new MergingIterator<Integer>(one_to_ten.iterator());
        x.add(evens.iterator());
        List<Integer> expected = Arrays.asList(1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10);

        testExpected( x, expected );
        Assert.assertFalse(x.hasNext());
    }

    @Test
    public void testMerging3() {
        logger.warn("Executing testMerging3");

        x = new MergingIterator<Integer>(Arrays.asList(bigs.iterator(), evens.iterator(), three_four_five.iterator()));
        List<Integer> expected = Arrays.asList(2, 3, 4, 4, 5, 6, 8, 10, 100, 200, 500);

        testExpected( x, expected );
        Assert.assertFalse(x.hasNext());
    }

    @Test
    public void testCollectLTE() {
        logger.warn("Executing testCollectLTE");

        x = new MergingIterator<Integer>(Arrays.asList(bigs.iterator(), evens.iterator(), three_four_five.iterator()));

        //for ( int i : x ) {
        //    System.out.printf("i=%d%n", i);
        //}

        logger.warn("Testing next element...");
        Assert.assertEquals( 2, (int)x.next());
        logger.warn("Testing allElementsLTE(2) = 2 ...");
        assertListsEqual( x.allElementsLTE(2), Arrays.asList(2) );

        logger.warn("Testing allElementsLTE(3) = 3 ...");
        Assert.assertEquals(3, (int)x.next());
        assertListsEqual( x.allElementsLTE(3), Arrays.asList(3) );

        logger.warn("Testing allElementsLTE(4) = 4, 4 ...");
        Assert.assertEquals(4, (int)x.next());
        assertListsEqual( x.allElementsLTE(4), Arrays.asList(4, 4) );

        logger.warn("Testing allElementsLTE(8) = 5, 6, 8 ...");
        assertListsEqual( x.allElementsLTE(8, false), Arrays.asList(5, 6, 8) );

        logger.warn("Testing remaining...");
        testExpected( x, Arrays.asList(10, 100, 200, 500) );
        Assert.assertFalse(x.hasNext());
    }

    private void testExpected( MergingIterator<Integer> mi, List<Integer> expected ) {
        for ( int expect : expected ) {
            int elt = mi.next();
            logger.warn(String.format("Testing elt=%d == expected=%d", elt, expect));
            Assert.assertEquals(expect, elt);
        }
    }

    private void assertListsEqual( Collection<Integer> l, List<Integer> expected ) {
        logger.warn("Size is " + l.size());
        Iterator<Integer> lit = l.iterator();

        for ( int expect : expected ) {
            Assert.assertTrue(lit.hasNext());
            int elt = lit.next();
            logger.warn(String.format("Testing elt=%d == expected=%d", elt, expect));
            Assert.assertEquals(expect, elt);
        }
    }
}