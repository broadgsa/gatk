/*
* Copyright 2012-2016 Broad Institute, Inc.
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

import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: 3/5/11
 * Time: 2:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class MWUnitTest extends BaseTest {
    private static double DELTA_PRECISION = 0.00001;

    @BeforeClass
    public void init() { }

    private static final MannWhitneyU rst = new MannWhitneyU();

    @DataProvider(name="rankSumTestData")
    public Object[][] dataProvider() {
        return new Object[][] {
                new Object[] {"test1", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,21}, 10d},
                new Object[] {"test2", new double[] {20,20,20,20,21}, new double[] {20,20,20,20,21}, 12.5d},
                new Object[] {"test3", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21}, 0.5d},
                new Object[] {"test4", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21,21,21,22,23,27}, 0.5d},
                new Object[] {"test5", new double[] {13,14,15,15,16,18,20,22,25,24,25,26,27,28,22,23,19,30,28,22,17},
                        new double[] {16,20,20,21,21,21,21,22,23,27,26,28,29,32,31,22,21,19,16,24,29},
                        180.5d},
                new Object[] {"test6", new double[] {13,14,15,15,16,18,13,14,15,15,16,18,13,14,15,15,16,18,13,14,15,15,16,18},
                        new double[] {21,22,23,27,26,28,29,21,22,23,27,26,28,29,21,22,23,27,26,28,29,21,22,23,27,26,28,29},
                        0d},
                new Object[] {"test7", new double[] {11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20},
                        new double[] {12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21},
                        145d},
                new Object[] {"test7", new double[] {11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20},
                        new double[] {12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21},
                        162d},
                new Object[] {"test8", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20}, 25d},
        };
    }

    @Test(dataProvider = "rankSumTestData")
    public void testSimpleU(String name, double[] series1, double[] series2, double U) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.TWO_SIDED);
        System.out.println("==========================================================");
        System.out.println("Series 1: " + Arrays.toString(series1));
        System.out.println("Series 2: " + Arrays.toString(series2));
        System.out.println("       U: " + test.getU());
        System.out.println("       Z: " + test.getZ());
        System.out.println("2-side p: " + test.getP());
        Assert.assertEquals(test.getU(), U, name);
    }

    @DataProvider(name="oneSidedPTestData")
    public Object[][] oneSidedDataProvider() {
        return new Object[][] {
                new Object[] {"test0", new double[] {0,0}, new double[] {1,1}, 0.083333333},
                new Object[] {"test1", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,21}, .25},
                new Object[] {"test2", new double[] {20,20,20,20,21}, new double[] {20,20,20,20,21}, .5},
                new Object[] {"test3", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21}, 0.00396825},
                new Object[] {"test4", new double[] {13,14,15,15,16}, new double[] {16,20,20,21,21,21,21,22,23,27}, 0.001469192},
                new Object[] {"test5", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20}, .5},
                new Object[] {"test6", new double[] {1,2,3,4,5}, new double[] {6,7,8,9,10}, 0.001984},
                new Object[] {"test7", new double[] {6,7,8,9,10}, new double[] {1,2,3,4,5}, 0.99801587},
                new Object[] {"test8", new double[] {16,20,20,21,21,21,21,22,23,27,16,20,20,21,21,21,21,22,23,27},
                        new double[] {22,23,27,16,20,20,21,21,21,21,22,23,27,60,60}, .08303102},
                new Object[] {"test9", new double[] {16,20,20,21,21,21,21,21,20},
                        new double[] {22,23,27,16,20,20,20,20,21}, 0.388204},
        };
    }

    @Test(dataProvider = "oneSidedPTestData")
    public void testOnesidedP(String name, double[] series1, double[] series2, double P) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.FIRST_DOMINATES);
        System.out.println("==========================================================");
        System.out.println("Series 1: " + Arrays.toString(series1));
        System.out.println("Series 2: " + Arrays.toString(series2));
        System.out.println("       U: " + test.getU());
        System.out.println("       Z: " + test.getZ());
        System.out.println("1-side p: " + test.getP());
        Assert.assertEquals(test.getP(), P, DELTA_PRECISION, name);
    }

    @DataProvider(name="oneSidedZTestData")
    public Object[][] oneSidedZDataProvider() {
        return new Object[][] {
                new Object[] {"test1", new double[] {20}, new double[] {20,20,20}, 0},
                new Object[] {"test2", new double[] {1}, new double[] {1,2,3}, -0.67448975},
                new Object[] {"test3", new double[] {1,2,3}, new double[] {3}, -0.67448975},
                new Object[] {"test4", new double[] {1,2,3}, new double[] {1}, 0.67448975},
                new Object[] {"test5", new double[] {3,3}, new double[] {1,2,3}, 1.036433},
                new Object[] {"test6", new double[] {20,20,20}, new double[] {20,20,20}, 0},
                new Object[] {"test7", new double[] {20,20,20,20,20}, new double[] {20,20,20,20,20,20,20,20,20,20,20,20,20}, 0},
                new Object[] {"test8", new double[] {1}, new double[] {2}, -0.67448975},
                new Object[] {"test9", new double[] {1}, new double[] {1}, 0},
                new Object[] {"test10", new double[] {60,70,70,60,60,60,60,60}, new double[] {60,60,60,60,60}, .91732119}
        };
    }

    @Test(dataProvider = "oneSidedZTestData")
    public void testOnesidedZ(String name, double[] series1, double[] series2, double Z) {
        MannWhitneyU.Result test = rst.test(series1, series2, MannWhitneyU.TestType.FIRST_DOMINATES);
        System.out.println("==========================================================");
        System.out.println("Series 1: " + Arrays.toString(series1));
        System.out.println("Series 2: " + Arrays.toString(series2));
        System.out.println("       U: " + test.getU());
        System.out.println("       Z: " + test.getZ());
        System.out.println("1-side p: " + test.getP());
        Assert.assertEquals(test.getZ(), Z, DELTA_PRECISION, name);
    }
}
