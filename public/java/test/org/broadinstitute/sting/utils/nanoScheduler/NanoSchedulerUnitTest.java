package org.broadinstitute.sting.utils.nanoScheduler;

import org.apache.log4j.BasicConfigurator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.SimpleTimer;
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
    private final static boolean debug = false;
    public static final int NANO_SCHEDULE_MAX_RUNTIME = 60000;

    private static class Map2x implements NSMapFunction<Integer, Integer> {
        @Override public Integer apply(Integer input) { return input * 2; }
    }

    private static class Map2xWithDelays extends Map2x {
        @Override public Integer apply(Integer input) {
            try {
                if ( input % 7 == 0 ) {
                    final int milliToSleep = (input % 10);
                    //System.out.printf("Sleeping %d millseconds%n", milliToSleep);
                    Thread.sleep(milliToSleep);
                }

                return input * 2;
            } catch ( InterruptedException ex ) {
                throw new RuntimeException(ex);
            }
        }
    }

    private static class ReduceSum implements NSReduceFunction<Integer, Integer> {
        int prevOne = Integer.MIN_VALUE;

        @Override public Integer apply(Integer one, Integer sum) {
            Assert.assertTrue(prevOne < one, "Reduce came in out of order.  Prev " + prevOne + " cur " + one);
            return one + sum;
        }
    }

    private static class ProgressCallback implements NSProgressFunction<Integer> {
        int callBacks = 0;

        @Override
        public void progress(Integer lastMapInput) {
            callBacks++;
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
        final boolean addDelays;

        public NanoSchedulerBasicTest(final int bufferSize, final int nThreads, final int start, final int end, final boolean addDelays) {
            super(NanoSchedulerBasicTest.class);
            this.bufferSize = bufferSize;
            this.nThreads = nThreads;
            this.start = start;
            this.end = end;
            this.expectedResult = sum2x(start, end);
            this.addDelays = addDelays;
            setName(String.format("%s nt=%d buf=%d start=%d end=%d sum=%d delays=%b",
                    getClass().getSimpleName(), nThreads, bufferSize, start, end, expectedResult, addDelays));
        }

        public Iterator<Integer> makeReader() {
            final List<Integer> ints = new ArrayList<Integer>();
            for ( int i = start; i < end; i++ )
                ints.add(i);
            return ints.iterator();
        }

        public int nExpectedCallbacks() {
            int nElements = Math.max(end - start, 0);
            return nElements / bufferSize;
        }

        public Map2x makeMap() { return addDelays ? new Map2xWithDelays() : new Map2x(); }
        public Integer initReduce() { return 0; }
        public ReduceSum makeReduce() { return new ReduceSum(); }

        public NanoScheduler<Integer, Integer, Integer> makeScheduler() {
            final NanoScheduler <Integer, Integer, Integer> nano;
            if ( bufferSize == -1 )
                nano = new NanoScheduler<Integer, Integer, Integer>(nThreads);
            else
                nano = new NanoScheduler<Integer, Integer, Integer>(bufferSize, nThreads);

            nano.setDebug(debug);
            return nano;
        }
    }

    static NanoSchedulerBasicTest exampleTest = null;
    @DataProvider(name = "NanoSchedulerBasicTest")
    public Object[][] createNanoSchedulerBasicTest() {
//        for ( final int bufferSize : Arrays.asList(1, 10) ) {
//            for ( final int nt : Arrays.asList(1, 2, 4) ) {
//                for ( final int start : Arrays.asList(0) ) {
//                    for ( final int end : Arrays.asList(0, 1, 2) ) {
//                        exampleTest = new NanoSchedulerBasicTest(bufferSize, nt, start, end, false);
//                    }
//                }
//            }
//        }

        for ( final int bufferSize : Arrays.asList(-1, 1, 10, 100) ) {
            for ( final int nt : Arrays.asList(1, 2, 4) ) {
                for ( final int start : Arrays.asList(0) ) {
                    for ( final int end : Arrays.asList(0, 1, 2, 11, 100, 10000, 100000) ) {
                        for ( final boolean addDelays : Arrays.asList(true, false) ) {
                            if ( end < 1000 )
                                exampleTest = new NanoSchedulerBasicTest(bufferSize, nt, start, end, addDelays);
                        }
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
        final SimpleTimer timer = new SimpleTimer().start();
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = test.makeScheduler();

        final ProgressCallback callback = new ProgressCallback();
        nanoScheduler.setProgressFunction(callback);

        if ( test.bufferSize > -1 )
            Assert.assertEquals(nanoScheduler.getBufferSize(), test.bufferSize, "bufferSize argument");
        Assert.assertEquals(nanoScheduler.getnThreads(), test.nThreads, "nThreads argument");

        final Integer sum = nanoScheduler.execute(test.makeReader(), test.makeMap(), test.initReduce(), test.makeReduce());
        Assert.assertNotNull(sum);
        Assert.assertEquals((int)sum, test.expectedResult, "NanoScheduler sum not the same as calculated directly");

        Assert.assertTrue(callback.callBacks >= test.nExpectedCallbacks(), "Not enough callbacks detected.  Expected at least " + test.nExpectedCallbacks() + " but saw only " + callback.callBacks);
        nanoScheduler.shutdown();

        // TODO -- need to enable only in the case where there's serious time spend in
        // TODO -- read /map / reduce, otherwise the "outside" timer doesn't add up
        final double myTimeEstimate = timer.getElapsedTime();
        final double tolerance = 0.1;
        if ( false && myTimeEstimate > 0.1 ) {
            Assert.assertTrue(nanoScheduler.getTotalRuntime() > myTimeEstimate * tolerance,
                    "NanoScheduler said that the total runtime was " + nanoScheduler.getTotalRuntime()
                            + " but the overall test time was " + myTimeEstimate + ", beyond our tolerance factor of "
                            + tolerance);
        }
    }

    @Test(enabled = true, dataProvider = "NanoSchedulerBasicTest", dependsOnMethods = "testMultiThreadedNanoScheduler", timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testNanoSchedulerInLoop(final NanoSchedulerBasicTest test) throws InterruptedException {
        if ( test.bufferSize > 1) {
            logger.warn("Running " + test);

            final NanoScheduler<Integer, Integer, Integer> nanoScheduler = test.makeScheduler();

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
        org.apache.log4j.Logger logger = org.apache.log4j.Logger.getRootLogger();
        BasicConfigurator.configure();
        logger.setLevel(org.apache.log4j.Level.DEBUG);

        final NanoSchedulerBasicTest test = new NanoSchedulerBasicTest(1000, Integer.valueOf(args[0]), 0, Integer.valueOf(args[1]), false);
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler =
                new NanoScheduler<Integer, Integer, Integer>(test.bufferSize, test.nThreads);
        nanoScheduler.setDebug(true);

        final Integer sum = nanoScheduler.execute(test.makeReader(), test.makeMap(), test.initReduce(), test.makeReduce());
        System.out.printf("Sum = %d, expected =%d%n", sum, test.expectedResult);
        nanoScheduler.shutdown();
    }
}
