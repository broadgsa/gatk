package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.utils.SimpleTimer;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Iterator;
import java.util.concurrent.BlockingQueue;

/**
 * Producer Thread that reads input values from an inputReads and puts them into a BlockingQueue
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

    public InputProducer(final Iterator<InputType> inputReader,
                         final SimpleTimer inputTimer,
                         final BlockingQueue<InputValue> outputQueue) {
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( outputQueue == null ) throw new IllegalArgumentException("OutputQueue cannot be null");

        this.inputReader = inputReader;
        this.inputTimer = inputTimer;
        this.outputQueue = outputQueue;
    }

    public void run() {
        try {
            while ( inputReader.hasNext() ) {
                if ( inputTimer != null ) inputTimer.restart();
                final InputType input = inputReader.next();
                if ( inputTimer != null ) inputTimer.stop();
                outputQueue.put(new InputValue(input));
            }

            // add the EOF object so our consumer knows we are done in all inputs
            outputQueue.put(new InputValue());
        } catch (InterruptedException ex) {
            throw new ReviewedStingException("got execution exception", ex);
        }
    }

    /**
     * Helper class that contains a read value suitable for EOF marking in a BlockingQueue
     */
    class InputValue extends BlockingQueueValue<InputType> {
        private InputValue(InputType datum) { super(datum); }
        private InputValue() { }
    }
}
