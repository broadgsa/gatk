package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingDeque;

/**
 * UnitTests for the InputProducer
 *
 * User: depristo
 * Date: 8/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class InputProducerUnitTest extends BaseTest {
    @DataProvider(name = "InputProducerTest")
    public Object[][] createInputProducerTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int nElements : Arrays.asList(0, 1, 10, 100, 1000, 10000, 100000) ) {
            for ( final int queueSize : Arrays.asList(1, 10, 100) ) {
                tests.add(new Object[]{ nElements, queueSize });
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "InputProducerTest", timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME)
    public void testInputProducer(final int nElements, final int queueSize) throws InterruptedException {
        final List<Integer> elements = new ArrayList<Integer>(nElements);
        for ( int i = 0; i < nElements; i++ ) elements.add(i);

        final LinkedBlockingDeque<InputProducer<Integer>.InputValue> readQueue =
                new LinkedBlockingDeque<InputProducer<Integer>.InputValue>(queueSize);

        final InputProducer<Integer> ip = new InputProducer<Integer>(elements.iterator(), null, readQueue);

        final ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(ip);

        int lastValue = -1;
        int nRead = 0;
        while ( true ) {
            final int observedQueueSize = readQueue.size();
            Assert.assertTrue(observedQueueSize <= queueSize,
                    "Reader is enqueuing more elements " + observedQueueSize + " than allowed " + queueSize);

            final InputProducer<Integer>.InputValue value = readQueue.take();
            if ( value.isLast() ) {
                Assert.assertEquals(nRead, nElements, "Number of input values " + nRead + " not all that are expected " + nElements);
                Assert.assertEquals(readQueue.size(), 0, "Last queue element found but queue contains more values!");
                break;
            } else {
                Assert.assertTrue(lastValue < value.getValue(), "Read values coming out of order!");
                final int expected = lastValue + 1;
                Assert.assertEquals((int)value.getValue(), expected, "Value observed " + value.getValue() + " not equal to the expected value " + expected);
                nRead++;
                lastValue = value.getValue();
            }
        }
    }
}
