/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.nanoScheduler;

import org.apache.log4j.BasicConfigurator;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.SimpleTimer;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
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
    private final static boolean DEBUG = false;
    private final static boolean debug = false;
    public static final int NANO_SCHEDULE_MAX_RUNTIME = 30000;
    public static final int EXCEPTION_THROWING_TEST_TIMEOUT = 10000;

    private static class Map2x implements NSMapFunction<Integer, Integer> {
        @Override public Integer apply(Integer input) { return input * 2; }
    }

    private static void maybeDelayMe(final int input) {
        try {
            if ( input % 7 == 0 ) {
                final int milliToSleep = (input % 10);
                //System.out.printf("Sleeping %d millseconds%n", milliToSleep);
                Thread.sleep(milliToSleep);
            }
        } catch ( InterruptedException ex ) {
            throw new RuntimeException(ex);
        }
    }

    private static class Map2xWithDelays extends Map2x {
        @Override public Integer apply(Integer input) {
            maybeDelayMe(input);
            return input * 2;
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
            return nElements / bufferSize / NanoScheduler.UPDATE_PROGRESS_FREQ;
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
    static NanoSchedulerBasicTest exampleTestWithDelays = null;

    @BeforeSuite
    public void setUp() throws Exception {
        exampleTest = new NanoSchedulerBasicTest(10, 2, 1, 10, false);
        exampleTestWithDelays = new NanoSchedulerBasicTest(10, 2, 1, 10, true);
    }

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
                                new NanoSchedulerBasicTest(bufferSize, nt, start, end, addDelays);
                        }
                    }
                }
            }
        }

        return NanoSchedulerBasicTest.getTests(NanoSchedulerBasicTest.class);
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "NanoSchedulerBasicTest", timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testSingleThreadedNanoScheduler(final NanoSchedulerBasicTest test) throws InterruptedException {
        logger.warn("Running " + test);
        if ( test.nThreads == 1 )
            testNanoScheduler(test);
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "NanoSchedulerBasicTest", timeOut = NANO_SCHEDULE_MAX_RUNTIME, dependsOnMethods = "testSingleThreadedNanoScheduler")
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
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "NanoSchedulerBasicTest", dependsOnMethods = "testMultiThreadedNanoScheduler", timeOut = 2 * NANO_SCHEDULE_MAX_RUNTIME)
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

    @Test(enabled = true && ! DEBUG, timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testShutdown() throws InterruptedException {
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = new NanoScheduler<Integer, Integer, Integer>(1, 2);
        Assert.assertFalse(nanoScheduler.isShutdown(), "scheduler should be alive");
        nanoScheduler.shutdown();
        Assert.assertTrue(nanoScheduler.isShutdown(), "scheduler should be dead");
    }

    @Test(enabled = true && ! DEBUG, expectedExceptions = IllegalStateException.class, timeOut = NANO_SCHEDULE_MAX_RUNTIME)
    public void testShutdownExecuteFailure() throws InterruptedException {
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = new NanoScheduler<Integer, Integer, Integer>(1, 2);
        nanoScheduler.shutdown();
        nanoScheduler.execute(exampleTest.makeReader(), exampleTest.makeMap(), exampleTest.initReduce(), exampleTest.makeReduce());
    }

    @DataProvider(name = "NanoSchedulerInputExceptionTest")
    public Object[][] createNanoSchedulerInputExceptionTest() {
        List<Object[]> tests = new ArrayList<Object[]>();


        for ( final int bufSize : Arrays.asList(100) ) {
            for ( final int nThreads : Arrays.asList(8) ) {
                for ( final boolean addDelays : Arrays.asList(true, false) ) {
                    final NanoSchedulerBasicTest test = new NanoSchedulerBasicTest(bufSize, nThreads, 1, 1000000, false);
                    final int maxN = addDelays ? 1000 : 10000;
                    for ( int nElementsBeforeError = 0; nElementsBeforeError < maxN; nElementsBeforeError += Math.max(nElementsBeforeError / 10, 1) ) {
                        tests.add(new Object[]{nElementsBeforeError, test, addDelays});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, expectedExceptions = NullPointerException.class, timeOut = EXCEPTION_THROWING_TEST_TIMEOUT)
    public void testInputErrorIsThrown_NPE() throws InterruptedException {
        executeTestErrorThrowingInput(10, new NullPointerException(), exampleTest, false);
    }

    @Test(enabled = true, expectedExceptions = ReviewedGATKException.class, timeOut = EXCEPTION_THROWING_TEST_TIMEOUT)
    public void testInputErrorIsThrown_RSE() throws InterruptedException {
        executeTestErrorThrowingInput(10, new ReviewedGATKException("test"), exampleTest, false);
    }

    @Test(enabled = true, expectedExceptions = NullPointerException.class, dataProvider = "NanoSchedulerInputExceptionTest", timeOut = EXCEPTION_THROWING_TEST_TIMEOUT, invocationCount = 1)
    public void testInputRuntimeExceptionDoesntDeadlock(final int nElementsBeforeError, final NanoSchedulerBasicTest test, final boolean addDelays ) throws InterruptedException {
        executeTestErrorThrowingInput(nElementsBeforeError, new NullPointerException(), test, addDelays);
    }

    @Test(enabled = true, expectedExceptions = ReviewedGATKException.class, dataProvider = "NanoSchedulerInputExceptionTest", timeOut = EXCEPTION_THROWING_TEST_TIMEOUT, invocationCount = 1)
    public void testInputErrorDoesntDeadlock(final int nElementsBeforeError, final NanoSchedulerBasicTest test, final boolean addDelays ) throws InterruptedException {
        executeTestErrorThrowingInput(nElementsBeforeError, new Error(), test, addDelays);
    }

    private void executeTestErrorThrowingInput(final int nElementsBeforeError, final Throwable ex, final NanoSchedulerBasicTest test, final boolean addDelays) {
        logger.warn("executeTestErrorThrowingInput " + nElementsBeforeError + " ex=" + ex + " test=" + test + " addInputDelays=" + addDelays);
        final NanoScheduler<Integer, Integer, Integer> nanoScheduler = test.makeScheduler();
        nanoScheduler.execute(new ErrorThrowingIterator(nElementsBeforeError, ex, addDelays), test.makeMap(), test.initReduce(), test.makeReduce());
    }

    private static class ErrorThrowingIterator implements Iterator<Integer> {
        final int nElementsBeforeError;
        final boolean addDelays;
        int i = 0;
        final Throwable ex;

        private ErrorThrowingIterator(final int nElementsBeforeError, Throwable ex, boolean addDelays) {
            this.nElementsBeforeError = nElementsBeforeError;
            this.ex = ex;
            this.addDelays = addDelays;
        }

        @Override public boolean hasNext() { return true; }
        @Override public Integer next() {
            if ( i++ > nElementsBeforeError ) {
                if ( ex instanceof Error )
                    throw (Error)ex;
                else if ( ex instanceof RuntimeException )
                    throw (RuntimeException)ex;
                else
                    throw new RuntimeException("Bad exception " + ex);
            } else if ( addDelays ) {
                maybeDelayMe(i);
                return i;
            } else {
                return i;
            }
        }
        @Override public void remove() { throw new UnsupportedOperationException("x"); }
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
