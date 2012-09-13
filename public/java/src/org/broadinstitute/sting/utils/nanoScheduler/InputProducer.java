package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Producer Thread that reads input values from an inputReads and puts them into an output queue
 */
class InputProducer<InputType> implements Runnable {
    /**
     * The iterator we are using to get data from
     */
    final Iterator<InputType> inputReader;

    /**
     * Our timer (may be null) that we use to track our input costs
     */
    final SimpleTimer inputTimer;

    /**
     * Where we put our input values for consumption
     */
    final BlockingQueue<InputValue> outputQueue;

    /**
     * Have we read the last value from inputReader?
     *
     * Must be a local variable, as inputReader.hasNext() can actually end up doing a lot
     * of work, and the method getNumInputValues() is supposed to be called not in the
     * thread executing the reading of values but in the thread enqueuing results
     */
    boolean readLastValue = false;

    int nRead = 0;

    /**
     * A latch used to block threads that want to start up only when all of the values
     * in inputReader have been read by the thread executing run()
     */
    final CountDownLatch latch = new CountDownLatch(1);

    public InputProducer(final Iterator<InputType> inputReader,
                         final SimpleTimer inputTimer,
                         final BlockingQueue<InputValue> outputQueue) {
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( inputTimer == null ) throw new IllegalArgumentException("inputTimer cannot be null");
        if ( outputQueue == null ) throw new IllegalArgumentException("OutputQueue cannot be null");

        this.inputReader = inputReader;
        this.inputTimer = inputTimer;
        this.outputQueue = outputQueue;
    }

    /**
     * Returns the number of elements in the input stream, AFTER we've read all of the values.
     * If we haven't read them all yet, returns -1
     *
     * @return the total number of elements in input stream, or -1 if some are still to be read
     */
    public synchronized int getNumInputValues() {
        return allInputsHaveBeenRead() ? nRead : -1;
    }

    /**
     * Returns true if all of the elements have been read from the input stream
     *
     * @return true if all of the elements have been read from the input stream
     */
    public synchronized boolean allInputsHaveBeenRead() {
        return readLastValue;
    }

    /**
     * Read the next item from the input stream, if possible
     *
     * If the inputReader has values, returns them, otherwise return null.
     *
     * This method is synchronized, as it manipulates local state accessed across multiple threads.
     *
     * @return the next input stream value, or null if the stream contains no more elements
     * @throws InterruptedException
     */
    private synchronized InputType readNextItem() throws InterruptedException {
        inputTimer.restart();
        if ( ! inputReader.hasNext() ) {
            // we are done, mark ourselves as such and return null
            readLastValue = true;
            inputTimer.stop();
            return null;
        } else {
            // get the next value, and return it
            final InputType input = inputReader.next();
            inputTimer.stop();
            nRead++;
            return input;
        }
    }

    /**
     * Run this input producer, looping over all items in the input reader and
     * enqueueing them as InputValues into the outputQueue.  After the
     * end of the stream has been encountered, any threads waiting because
     * they called waitForDone() will be freed.
     */
    public void run() {
        try {
            while ( true ) {
                final InputValue inputValue = runOne();
                outputQueue.put(inputValue);
                if ( inputValue.isEOFMarker() )
                    break;
            }

            latch.countDown();
        } catch (InterruptedException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        }
    }

    protected InputValue runOne() throws InterruptedException {
        final InputType value = readNextItem();
        if ( value == null ) {
            // add the EOF object so our consumer knows we are done in all inputs
            return new InputValue();
        } else {
            // add the actual value
            return new InputValue(value);
        }
    }

    /**
     * Block until all of the items have been read from inputReader.
     *
     * Note that this call doesn't actually read anything.   You have to submit a thread
     * to actually execute run() directly.
     *
     * @throws InterruptedException
     */
    public void waitForDone() throws InterruptedException {
        latch.await();
    }

    /**
     * Helper class that contains a read value suitable for EOF marking in a BlockingQueue
     */
    class InputValue extends EOFMarkedValue<InputType> {
        private InputValue(InputType datum) { super(datum); }
        private InputValue() { }
    }
}
