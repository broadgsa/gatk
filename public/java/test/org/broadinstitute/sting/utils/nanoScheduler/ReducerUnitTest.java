package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MultiThreadedErrorTracker;
import org.broadinstitute.sting.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.TimeUnit;

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
            for ( final boolean setJobIDAtStart : Arrays.asList(true, false) ) {
                for ( final int nElements : Arrays.asList(0, 1, 3, 5) ) {
                    if ( groupSize < nElements ) {
                        for ( final List<MapResult<Integer>> jobs : Utils.makePermutations(makeJobs(nElements), nElements, false) ) {
                            tests.add(new Object[]{ new ListOfJobs(jobs), setJobIDAtStart, groupSize });
                        }
                    }
                }

                for ( final int nElements : Arrays.asList(10, 100, 1000, 10000, 100000, 1000000) ) {
                    if ( groupSize < nElements ) {
                        tests.add(new Object[]{ new ListOfJobs(makeJobs(nElements)), setJobIDAtStart, groupSize });
                    }
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
    public void testReducerThread(final List<MapResult<Integer>> jobs, final boolean setJobIDAtStart, final int groupSize) throws Exception {
        runTests(jobs, setJobIDAtStart, groupSize);
    }

    private void runTests( final List<MapResult<Integer>> allJobs, boolean setJobIDAtStart, int groupSize ) throws Exception {
        if ( groupSize == -1 )
            groupSize = allJobs.size();

        final PriorityBlockingQueue<MapResult<Integer>> mapResultsQueue = new PriorityBlockingQueue<MapResult<Integer>>();

        final List<List<MapResult<Integer>>> jobGroups = Utils.groupList(allJobs, groupSize);
        final ReduceSumTest reduce = new ReduceSumTest();
        final Reducer<Integer, Integer> reducer = new Reducer<Integer, Integer>(reduce, new MultiThreadedErrorTracker(), 0);

        final TestWaitingForFinalReduce waitingThread = new TestWaitingForFinalReduce(reducer, expectedSum(allJobs));
        final ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(waitingThread);

        int nJobsSubmitted = 0;
        int jobGroupCount = 0;
        final int lastJobGroupCount = jobGroups.size() - 1;
        setJobIDAtStart = setJobIDAtStart && groupSize == 1;

        for ( final List<MapResult<Integer>> jobs : jobGroups ) {
            //logger.warn("Processing job group " + jobGroupCount + " with " + jobs.size() + " jobs");
            for ( final MapResult<Integer> job : jobs ) {
                mapResultsQueue.add(job);
                nJobsSubmitted++;
            }

            if ( jobGroupCount == lastJobGroupCount ) {
                mapResultsQueue.add(new MapResult<Integer>());
                nJobsSubmitted++;
            }

            Assert.assertFalse(reducer.latchIsReleased(), "Latch should be closed at the start");

            if ( jobGroupCount == 0 && setJobIDAtStart ) {
                // only can do the setJobID if jobs cannot be submitted out of order
                reducer.setTotalJobCount(allJobs.size());
                Assert.assertFalse(reducer.latchIsReleased(), "Latch should be closed even after setting last job if we haven't processed anything");
            }

            final int nReduced = reducer.reduceAsMuchAsPossible(mapResultsQueue);
            Assert.assertTrue(nReduced <= nJobsSubmitted, "Somehow reduced more jobs than submitted");

            if ( setJobIDAtStart ) {
                final boolean submittedLastJob = jobGroupCount == lastJobGroupCount;
                Assert.assertEquals(reducer.latchIsReleased(), submittedLastJob,
                        "When last job is set, latch should only be released if the last job has been submitted");
            } else {
                Assert.assertEquals(reducer.latchIsReleased(), false, "When last job isn't set, latch should never be release");
            }

            jobGroupCount++;
        }

        if ( setJobIDAtStart )
            Assert.assertTrue(reducer.latchIsReleased(), "Latch should be released after reducing with last job id being set");
        else {
            Assert.assertFalse(reducer.latchIsReleased(), "Latch should be closed after reducing without last job id being set");
            reducer.setTotalJobCount(allJobs.size());
            Assert.assertTrue(reducer.latchIsReleased(), "Latch should be released after reducing after setting last job id ");
        }

        Assert.assertEquals(reduce.nRead, allJobs.size(), "number of read values not all of the values in the reducer queue");
        es.shutdown();
        es.awaitTermination(1, TimeUnit.HOURS);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    private void runSettingJobIDTwice() throws Exception {
        final PriorityBlockingQueue<MapResult<Integer>> mapResultsQueue = new PriorityBlockingQueue<MapResult<Integer>>();

        final Reducer<Integer, Integer> reducer = new Reducer<Integer, Integer>(new ReduceSumTest(), new MultiThreadedErrorTracker(), 0);

        reducer.setTotalJobCount(10);
        reducer.setTotalJobCount(15);
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
            try {
                final int observedSum = reducer.waitForFinalReduce();
                Assert.assertEquals(observedSum, expectedSum, "Reduce didn't sum to expected value");
            } catch ( InterruptedException ex ) {
                Assert.fail("Got interrupted");
            }
        }
    }
}