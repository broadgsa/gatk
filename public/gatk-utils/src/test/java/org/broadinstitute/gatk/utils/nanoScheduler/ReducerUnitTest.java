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

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.MultiThreadedErrorTracker;
import org.broadinstitute.gatk.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.*;

/**
 * UnitTests for Reducer
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReducerUnitTest extends BaseTest {
    @DataProvider(name = "ReducerThreadTest")
    public Object[][] createReducerThreadTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int groupSize : Arrays.asList(-1, 1, 5, 50, 500, 5000, 50000) ) {
            for ( final int nElements : Arrays.asList(0, 1, 3, 5) ) {
                if ( groupSize < nElements ) {
                    for ( final List<MapResult<Integer>> jobs : Utils.makePermutations(makeJobs(nElements), nElements, false) ) {
                        tests.add(new Object[]{ new ListOfJobs(jobs), groupSize });
                    }
                }
            }

            for ( final int nElements : Arrays.asList(10, 100, 1000, 10000, 100000, 1000000) ) {
                if ( groupSize < nElements ) {
                    tests.add(new Object[]{ new ListOfJobs(makeJobs(nElements)), groupSize });
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private static class ListOfJobs extends ArrayList<MapResult<Integer>> {
        private ListOfJobs(Collection<? extends MapResult<Integer>> c) {
            super(c);
        }

        @Override
        public String toString() {
            if ( size() < 10 )
                return super.toString();
            else
                return "JobList of " + size();
        }
    }

    private static List<MapResult<Integer>> makeJobs(final int nElements) {
        List<MapResult<Integer>> jobs = new ArrayList<MapResult<Integer>>(nElements);
        for ( int i = 0; i < nElements; i++ ) {
            jobs.add(new MapResult<Integer>(i, i));
        }
        return jobs;
    }

    private int expectedSum(final List<MapResult<Integer>> jobs) {
        int sum = 0;
        for ( final MapResult<Integer> job : jobs )
            sum += job.getValue();
        return sum;
    }

    @Test(enabled = true, dataProvider = "ReducerThreadTest", timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME)
    public void testReducerThread(final List<MapResult<Integer>> allJobs, int groupSize) throws Exception {
        if ( groupSize == -1 )
            groupSize = allJobs.size();

        final MapResultsQueue<Integer> mapResultsQueue = new MapResultsQueue<Integer>();

        final List<List<MapResult<Integer>>> jobGroups = Utils.groupList(allJobs, groupSize);
        final ReduceSumTest reduce = new ReduceSumTest();
        final Reducer<Integer, Integer> reducer = new Reducer<Integer, Integer>(reduce, new MultiThreadedErrorTracker(), 0);

        final TestWaitingForFinalReduce waitingThread = new TestWaitingForFinalReduce(reducer, expectedSum(allJobs));
        final ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(waitingThread);

        int lastJobID = -1;
        int nJobsSubmitted = 0;
        int jobGroupCount = 0;
        final int lastJobGroupCount = jobGroups.size() - 1;

        for ( final List<MapResult<Integer>> jobs : jobGroups ) {
            //logger.warn("Processing job group " + jobGroupCount + " with " + jobs.size() + " jobs");
            for ( final MapResult<Integer> job : jobs ) {
                lastJobID = Math.max(lastJobID, job.getJobID());
                mapResultsQueue.put(job);
                nJobsSubmitted++;
            }

            if ( jobGroupCount == lastJobGroupCount ) {
                mapResultsQueue.put(new MapResult<Integer>(lastJobID+1));
                nJobsSubmitted++;
            }

            final int nReduced = reducer.reduceAsMuchAsPossible(mapResultsQueue, true);
            Assert.assertTrue(nReduced <= nJobsSubmitted, "Somehow reduced more jobs than submitted");

            jobGroupCount++;
        }

        Assert.assertEquals(reduce.nRead, allJobs.size(), "number of read values not all of the values in the reducer queue");
        es.shutdown();
        es.awaitTermination(1, TimeUnit.HOURS);
    }

    @Test(timeOut = 1000, invocationCount = 100)
    private void testNonBlockingReduce() throws Exception {
        final Reducer<Integer, Integer> reducer = new Reducer<Integer, Integer>(new ReduceSumTest(), new MultiThreadedErrorTracker(), 0);
        final MapResultsQueue<Integer> mapResultsQueue = new MapResultsQueue<Integer>();
        mapResultsQueue.put(new MapResult<Integer>(0, 0));
        mapResultsQueue.put(new MapResult<Integer>(1, 1));

        final CountDownLatch latch = new CountDownLatch(1);
        final ExecutorService es = Executors.newSingleThreadExecutor();

        es.submit(new Runnable() {
            @Override
            public void run() {
                reducer.acquireReduceLock(true);
                latch.countDown();
            }
        });

        latch.await();
        final int nReduced = reducer.reduceAsMuchAsPossible(mapResultsQueue, false);
        Assert.assertEquals(nReduced, 0, "The reducer lock was already held but we did some work");
        es.shutdown();
        es.awaitTermination(1, TimeUnit.HOURS);
    }

    @Test(timeOut = 10000, invocationCount = 100)
    private void testBlockingReduce() throws Exception {
        final Reducer<Integer, Integer> reducer = new Reducer<Integer, Integer>(new ReduceSumTest(), new MultiThreadedErrorTracker(), 0);
        final MapResultsQueue<Integer> mapResultsQueue = new MapResultsQueue<Integer>();
        mapResultsQueue.put(new MapResult<Integer>(0, 0));
        mapResultsQueue.put(new MapResult<Integer>(1, 1));

        final CountDownLatch latch = new CountDownLatch(1);
        final ExecutorService es = Executors.newSingleThreadExecutor();

        es.submit(new Runnable() {
            @Override
            public void run() {
                reducer.acquireReduceLock(true);
                latch.countDown();
                try {
                    Thread.sleep(100);
                } catch ( InterruptedException e ) {
                    ;
                } finally {
                    reducer.releaseReduceLock();
                }
            }
        });

        latch.await();
        final int nReduced = reducer.reduceAsMuchAsPossible(mapResultsQueue, true);
        Assert.assertEquals(nReduced, 2, "The reducer should have blocked until the lock was freed and reduced 2 values");
        es.shutdown();
        es.awaitTermination(1, TimeUnit.HOURS);
    }


    public class ReduceSumTest implements NSReduceFunction<Integer, Integer> {
        int nRead = 0;
        int lastValue = -1;

        @Override public Integer apply(Integer one, Integer sum) {
            Assert.assertTrue(lastValue < one, "Reduce came in out of order.  Prev " + lastValue + " cur " + one);

            Assert.assertTrue(lastValue < one, "Read values coming out of order!");
            final int expected = lastValue + 1;
            Assert.assertEquals((int)one, expected, "Value observed " + one + " not equal to the expected value " + expected);
            nRead++;
            lastValue = expected;

            return one + sum;
        }
    }

    final static class TestWaitingForFinalReduce implements Runnable {
        final Reducer<Integer, Integer> reducer;
        final int expectedSum;

        TestWaitingForFinalReduce(Reducer<Integer, Integer> reducer, final int expectedSum) {
            this.reducer = reducer;
            this.expectedSum = expectedSum;
        }

        @Override
        public void run() {
            final int observedSum = reducer.getReduceResult();
            Assert.assertEquals(observedSum, expectedSum, "Reduce didn't sum to expected value");
        }
    }
}