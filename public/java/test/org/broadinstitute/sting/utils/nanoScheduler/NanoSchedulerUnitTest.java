package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * UnitTests for the NanoScheduler
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class NanoSchedulerUnitTest extends BaseTest {
    private class Map2x implements MapFunction<Integer, Integer> {
        @Override public Integer apply(Integer input) { return input * 2; }
    }

    private class ReduceSum implements ReduceFunction<Integer, Integer> {
        @Override public Integer init() { return 0; }
        @Override public Integer apply(Integer one, Integer sum) { return one + sum; }
    }

    private static int sum2x(final int start, final int end) {
        int sum = 0;
        for ( int i = start; i < end; i++ )
            sum += 2 * i;
        return sum;
    }

    private class NanoSchedulerBasicTest extends TestDataProvider {
        final int bufferSize, nThreads, start, end, expectedResult;

        public NanoSchedulerBasicTest(final int bufferSize, final int nThreads, final int start, final int end) {
            super(NanoSchedulerBasicTest.class);
            this.bufferSize = bufferSize;
            this.nThreads = nThreads;
            this.start = start;
            this.end = end;
            this.expectedResult = sum2x(start, end);
            setName(String.format("%s nt=%d buf=%d start=%d end=%d sum=%d",
                    getClass().getSimpleName(), nThreads, bufferSize, start, end, expectedResult));
        }

        public Iterator<Integer> makeReader() {
            final List<Integer> ints = new ArrayList<Integer>();
            for ( int i = start; i < end; i++ )
                ints.add(i);
            return ints.iterator();
        }

        public Map2x makeMap() { return new Map2x(); }
        public ReduceSum makeReduce() { return new ReduceSum(); }
    }

    @DataProvider(name = "NanoSchedulerBasicTest")
    public Object[][] createNanoSchedulerBasicTest() {
        for ( final int bufferSize : Arrays.asList(1, 10, 10000, 1000000) ) {
            for ( final int nt : Arrays.asList(1, 2, 4, 8, 16, 32) ) {
                for ( final int start : Arrays.asList(0) ) {
                    for ( final int end : Arrays.asList(1, 2, 11, 1000000) ) {
                        new NanoSchedulerBasicTest(bufferSize,  nt, start, end);
                    }
                }
            }
        }

        return NanoSchedulerBasicTest.getTests(NanoSchedulerBasicTest.class);
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", timeOut = 2000)
    public void testNanoSchedulerBasicTest(final NanoSchedulerBasicTest test) throws InterruptedException {
        logger.warn("Running " + test);
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler =
                new NanoScheduler<Integer, Integer, Integer>(test.bufferSize, test.nThreads,
                        test.makeReader(), test.makeMap(), test.makeReduce());
        final Integer sum = nanoScheduler.execute();
        Assert.assertNotNull(sum);
        Assert.assertEquals((int)sum, test.expectedResult, "NanoScheduler sum not the same as calculated directly");
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", timeOut = 10000, dependsOnMethods = "testNanoSchedulerBasicTest")
    public void testNanoSchedulerInLoop(final NanoSchedulerBasicTest test) throws InterruptedException {
        logger.warn("Running " + test);
        for ( int i = 0; i < 10; i++ ) {
            testNanoSchedulerBasicTest(test);
        }
    }
}
