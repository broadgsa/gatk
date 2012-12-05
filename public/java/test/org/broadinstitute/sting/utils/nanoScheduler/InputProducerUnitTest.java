package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.MultiThreadedErrorTracker;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.Semaphore;

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

        final InputProducer<Integer> ip = new InputProducer<Integer>(elements.iterator(), new MultiThreadedErrorTracker(), readQueue);

        final ExecutorService es = Executors.newSingleThreadExecutor();

        Assert.assertFalse(ip.allInputsHaveBeenRead(), "InputProvider said that all inputs have been read, but I haven't started reading yet");
        Assert.assertEquals(ip.getNumInputValues(), -1, "InputProvider told me that the queue was done, but I haven't started reading yet");

        es.submit(ip);

        int lastValue = -1;
        int nRead = 0;
        while ( true ) {
            final int nTotalElements = ip.getNumInputValues();
            final int observedQueueSize = readQueue.size();
            Assert.assertTrue(observedQueueSize <= queueSize,
                    "Reader is enqueuing more elements " + observedQueueSize + " than allowed " + queueSize);

            if ( nRead + observedQueueSize < nElements )
                Assert.assertEquals(nTotalElements, -1, "getNumInputValues should have returned -1 with not all elements read");
            // note, cannot test else case because elements input could have emptied between calls

            final InputProducer<Integer>.InputValue value = readQueue.take();
            if ( value.isEOFMarker() ) {
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

        Assert.assertTrue(ip.allInputsHaveBeenRead(), "InputProvider said that all inputs haven't been read, but I read them all");
        Assert.assertEquals(ip.getNumInputValues(), nElements, "Wrong number of total elements getNumInputValues");
        es.shutdownNow();
    }

    @Test(enabled = true, dataProvider = "InputProducerTest", timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME)
    public void testInputProducerLocking(final int nElements, final int queueSize) throws InterruptedException {
        final List<Integer> elements = new ArrayList<Integer>(nElements);
        for ( int i = 0; i < nElements; i++ ) elements.add(i);

        final LinkedBlockingDeque<InputProducer<Integer>.InputValue> readQueue =
                new LinkedBlockingDeque<InputProducer<Integer>.InputValue>();

        final InputProducer<Integer> ip = new InputProducer<Integer>(elements.iterator(), new MultiThreadedErrorTracker(), readQueue);

        final ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(ip);

        ip.waitForDone();

        Assert.assertEquals(ip.getNumInputValues(), nElements, "InputProvider told me that the queue was done, but I haven't started reading yet");
        Assert.assertEquals(readQueue.size(), nElements + 1, "readQueue should have had all elements read into it");
    }

    final static class BlockingIterator<T> implements Iterator<T> {
        final Semaphore blockNext = new Semaphore(0);
        final Semaphore blockOnNext = new Semaphore(0);
        final Iterator<T> underlyingIterator;

        BlockingIterator(Iterator<T> underlyingIterator) {
            this.underlyingIterator = underlyingIterator;
        }

        public void allowNext() {
            blockNext.release(1);
        }

        public void blockTillNext() throws InterruptedException {
            blockOnNext.acquire(1);
        }

        @Override
        public boolean hasNext() {
            return underlyingIterator.hasNext();
        }

        @Override
        public T next() {
            try {
                blockNext.acquire(1);
                T value = underlyingIterator.next();
                blockOnNext.release(1);
                return value;
            } catch (InterruptedException ex) {
                throw new RuntimeException(ex);
            }
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("x");
        }
    }
}
