// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.

import org.junit.Assert;
import org.junit.Test;
import org.junit.Before;
import org.broadinstitute.sting.BaseTest;
import java.util.Arrays;

/**
 * Basic unit test for RecalData
 */
public class ExpandingArrayListTest extends BaseTest {
    ExpandingArrayList<Integer> empty, initCap10, hasOne, hasTen;

    @Before
    public void before() {
        empty = new ExpandingArrayList<Integer>();

        initCap10 = new ExpandingArrayList<Integer>(10);

        hasOne = new ExpandingArrayList<Integer>();
        hasOne.add(1);

        hasTen = new ExpandingArrayList<Integer>();
        hasTen.addAll(Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
    }

    @Test
    public void testBasicSizes() {
        logger.warn("Executing testBasicSizes");

        Assert.assertEquals(empty.size(), 0);
        Assert.assertEquals(initCap10.size(), 0);
        Assert.assertEquals(hasOne.size(), 1);
        Assert.assertEquals(hasTen.size(), 10);
    }

    @Test
    public void testTenElements() {
        logger.warn("Executing testTenElements");

        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals((int)hasTen.get(i), i+1);
        }
    }

    @Test
    public void testSettingTenElements() {
        logger.warn("Executing testSettingTenElements");

        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals((int)hasTen.set(i, 2*i), i+1);
        }

        Assert.assertEquals(hasTen.size(), 10);
        for ( int i = 0; i < 10; i++ ) {
            Assert.assertEquals((int)hasTen.get(i), 2*i);
        }
    }

    @Test
    public void testAdd() {
        logger.warn("Executing testAdd");

        Assert.assertEquals(empty.size(), 0);
        empty.add(1);
        Assert.assertEquals(empty.size(), 1);
        Assert.assertEquals((int)empty.get(0), 1);
        empty.add(2);
        Assert.assertEquals(empty.size(), 2);
        Assert.assertEquals((int)empty.get(1), 2);        
    }

    @Test
    public void testSet1() {
        logger.warn("Executing testSet1");

        Assert.assertEquals(empty.size(), 0);
        empty.set(0, 1);
        Assert.assertEquals(empty.size(), 1);
        Assert.assertEquals((int)empty.get(0), 1);

        empty.set(1, 2);
        Assert.assertEquals(empty.size(), 2);
        Assert.assertEquals((int)empty.get(1), 2);

        // doesn't expand
        empty.set(0, 3);
        Assert.assertEquals(empty.size(), 2);
        Assert.assertEquals((int)empty.get(0), 3);
    }

    @Test
    public void testSetExpanding() {
        logger.warn("Executing testSetExpanding");

        Assert.assertEquals(empty.size(), 0);
        empty.set(3, 1);
        Assert.assertEquals(empty.size(), 4);
        Assert.assertEquals(empty.get(0), null);
        Assert.assertEquals(empty.get(1), null);
        Assert.assertEquals(empty.get(2), null);
        Assert.assertEquals((int)empty.get(3), 1);
    }

    @Test
    public void testSetExpandingReset() {
        logger.warn("Executing testSetExpandingReset");

        Assert.assertEquals(empty.size(), 0);
        empty.set(3, 3);
        empty.set(2, 2);
        empty.set(1, 1);
        empty.set(0, 0);
        Assert.assertEquals(empty.size(), 4);
        for ( int i = 0; i < 4; i++ )
            Assert.assertEquals((int)empty.get(i), i);
    }

    @Test
    public void testSetExpandingBig() {
        logger.warn("Executing testSetExpandingBig");

        Assert.assertEquals(empty.size(), 0);
        empty.set(1000, 1000);
        Assert.assertEquals(empty.size(), 1001);
        for ( int i = 0; i < 1000; i++ )
            Assert.assertEquals(empty.get(i), null);
        Assert.assertEquals((int)empty.get(1000), 1000);
    }

    @Test (expected=IndexOutOfBoundsException.class )
    public void testSetBadGetNegative() {
        logger.warn("Executing testSetBadGetNegative");
        empty.get(-1);
    }

    @Test
    public void testSetBadGetPost() {
        logger.warn("Executing testSetBadGetPost");
        empty.set(1, 1);
        Assert.assertEquals(empty.get(0), null);
        Assert.assertEquals((int)empty.get(1), 1);
        Assert.assertEquals(empty.get(2), null);
    }
}