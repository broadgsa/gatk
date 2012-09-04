package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * UnitTests for the NanoScheduler
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class NanoSchedulerUnitTest extends BaseTest {
    public static final int NANO_SCHEDULE_MAX_RUNTIME = 60000;

    private static class Map2x implements NanoSchedulerMapFunction<Integer, Integer> {
        @Override public Integer apply(Integer input) { return input * 2; }
    }

    private static class ReduceSum implements NanoSchedulerReduceFunction<Integer, Integer> {
        int prevOne = Integer.MIN_VALUE;

        @Override public Integer apply(Integer one, Integer sum) {
            Assert.assertTrue(prevOne < one, "Reduce came in out of order.  Prev " + prevOne + " cur " + one);
            return one + sum;
        }
    }

    private static int sum2x(final int start, final int end) {
        int sum = 0;
        for ( int i = start; i < end; i++ )
            sum += 2 * i;
        return sum;
    }

    private static class NanoSchedulerBasicTest extends TestDataProvider {
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
        public Integer initReduce() { return 0; }
        public ReduceSum makeReduce() { return new ReduceSum(); }
    }

    static NanoSchedulerBasicTest exampleTest = null;
    @DataProvider(name = "NanoSchedulerBasicTest")
    public Object[][] createNanoSchedulerBasicTest() {
        for ( final int bufferSize : Arrays.asList(1, 10, 1000, 1000000) ) {
            for ( final int nt : Arrays.asList(1, 2, 4) ) {
                for ( final int start : Arrays.asList(0) ) {
                    for ( final int end : Arrays.asList(1, 2, 11, 10000, 100000) ) {
                        exampleTest = new NanoSchedulerBasicTest(bufferSize, nt, start, end);
                    }
                }
            }
        }

        return NanoSchedulerBasicTest.getTests(NanoSchedulerBasicTest.class);
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testSingleThreadedNanoScheduler(final NanoSchedulerBasicTest test) throws InterruptedException {
        logger.warn("Running " + test);
        if ( test.nThreads == 1 )
            testNanoScheduler(test);
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", timeOut = NANO_SCHEDULE_MAX_RUNTIME, dependsOnMethods = "testSingleThreadedNanoScheduler")
    public void testMultiThreadedNanoScheduler(final NanoSchedulerBasicTest test) throws InterruptedException {
        logger.warn("Running " + test);
        if ( test.nThreads >= 1 )
            testNanoScheduler(test);
    }

    private void testNanoScheduler(final NanoSchedulerBasicTest test) throws InterruptedException {
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler =
                new NanoScheduler<Integer, Integer, Integer>(test.bufferSize, test.nThreads);

        Assert.assertEquals(nanoScheduler.getBufferSize(), test.bufferSize, "bufferSize argument");
        Assert.assertEquals(nanoScheduler.getnThreads(), test.nThreads, "nThreads argument");

        final Integer sum = nanoScheduler.execute(test.makeReader(), test.makeMap(), test.initReduce(), test.makeReduce());
        Assert.assertNotNull(sum);
        Assert.assertEquals((int)sum, test.expectedResult, "NanoScheduler sum not the same as calculated directly");
        nanoScheduler.shutdown();
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", dependsOnMethods = "testMultiThreadedNanoScheduler", timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testNanoSchedulerInLoop(final NanoSchedulerBasicTest test) throws InterruptedException {
        if ( test.bufferSize > 1) {
            logger.warn("Running " + test);

            final NanoScheduler<Integer, Integer, Integer> nanoScheduler =
                    new NanoScheduler<Integer, Integer, Integer>(test.bufferSize, test.nThreads);

            // test reusing the scheduler
            for ( int i = 0; i < 10; i++ ) {
                final Integer sum = nanoScheduler.execute(test.makeReader(), test.makeMap(), test.initReduce(), test.makeReduce());
                Assert.assertNotNull(sum);
                Assert.assertEquals((int)sum, test.expectedResult, "NanoScheduler sum not the same as calculated directly");
            }

            nanoScheduler.shutdown();
        }
    }

    @Test(timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testShutdown() throws InterruptedException {
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = new NanoScheduler<Integer, Integer, Integer>(1, 2);
        Assert.assertFalse(nanoScheduler.isShutdown(), "scheduler should be alive");
        nanoScheduler.shutdown();
        Assert.assertTrue(nanoScheduler.isShutdown(), "scheduler should be dead");
    }

    @Test(expectedExceptions = IllegalStateException.class, timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testShutdownExecuteFailure() throws InterruptedException {
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = new NanoScheduler<Integer, Integer, Integer>(1, 2);
        nanoScheduler.shutdown();
        nanoScheduler.execute(exampleTest.makeReader(), exampleTest.makeMap(), exampleTest.initReduce(), exampleTest.makeReduce());
    }

    public static void main(String [ ] args) {
        final NanoSchedulerBasicTest test = new NanoSchedulerBasicTest(1000, Integer.valueOf(args[0]), 0, Integer.valueOf(args[1]));
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler =
                new NanoScheduler<Integer, Integer, Integer>(test.bufferSize, test.nThreads);

        final Integer sum = nanoScheduler.execute(test.makeReader(), test.makeMap(), test.initReduce(), test.makeReduce());
        System.out.printf("Sum = %d, expected =%d%n", sum, test.expectedResult);
    }
}
