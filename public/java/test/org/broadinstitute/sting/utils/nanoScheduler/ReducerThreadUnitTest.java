package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

/**
 * UnitTests for the InputProducer
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReducerThreadUnitTest extends BaseTest {
    @DataProvider(name = "ReducerThreadTest")
    public Object[][] createReducerThreadTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int nElements : Arrays.asList(0, 1, 10, 100, 1000, 10000, 100000) ) {
            tests.add(new Object[]{ nElements });
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "ReducerThreadTest", timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME)
    public void testReducerThreadTest(final int nElements) throws Exception {
        List<Integer> values = new ArrayList<Integer>(nElements);
        List<Integer> jobIDs = new ArrayList<Integer>(nElements);
        for ( int i = 0; i < nElements; i++ ) {
            values.add(i);
            jobIDs.add(i);
        }

        runTests(values, jobIDs);
    }

    @Test(enabled = true, timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME, expectedExceptions = ExecutionException.class)
    public void testReducerThreadTestByJobOrder() throws Exception {
        runTests(Arrays.asList(0, 1, 2), Arrays.asList(1, 3, 2));
    }

    private void runTests( final List<Integer> mapValues, final List<Integer> jobIDs) throws Exception {
        final LinkedBlockingDeque<Future<MapResult<Integer>>> mapResultsQueue =
                new LinkedBlockingDeque<Future<MapResult<Integer>>>(mapValues.size()+1);

        for ( int i = 0; i < mapValues.size(); i++ ) {
            final int value = mapValues.get(i);
            final int jobID = jobIDs.get(i);
            final MapResult<Integer> mapResult = new MapResult<Integer>(value, jobID);
            mapResultsQueue.add(new FutureValue<MapResult<Integer>>(mapResult));
        }
        mapResultsQueue.add(new FutureValue<MapResult<Integer>>(new MapResult<Integer>()));

        final ReduceSumTest reduce = new ReduceSumTest(mapResultsQueue);
        final ReducerThread<Integer, Integer> thread
                = new ReducerThread<Integer, Integer>(reduce, null, 0, mapResultsQueue);

        final ExecutorService es = Executors.newSingleThreadExecutor();
        final Future<Integer> value = es.submit(thread);
        value.get();

        Assert.assertEquals(reduce.nRead, mapValues.size());
    }

    public class ReduceSumTest implements NSReduceFunction<Integer, Integer> {
        final LinkedBlockingDeque<Future<MapResult<Integer>>> mapResultsQueue;
        int nRead = 0;
        int lastValue = -1;

        public ReduceSumTest(LinkedBlockingDeque<Future<MapResult<Integer>>> mapResultsQueue) {
            this.mapResultsQueue = mapResultsQueue;
        }

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
}
