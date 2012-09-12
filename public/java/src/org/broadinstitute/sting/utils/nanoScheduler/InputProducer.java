package org.broadinstitute.sting.utils.nanoScheduler;

import org.broadinstitute.sting.utils.SimpleTimer;

import java.util.Iterator;

/**
 * Producer Thread that reads input values from an inputReads and puts them into a BlockingQueue
 */
class InputProducer<InputType> {
    /**
     * The iterator we are using to get data from
     */
    final Iterator<InputType> inputReader;

    /**
     * Our timer (may be null) that we use to track our input costs
     */
    final SimpleTimer inputTimer;

    public InputProducer(final Iterator<InputType> inputReader,
                         final SimpleTimer inputTimer) {
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( inputTimer == null ) throw new IllegalArgumentException("inputTimer cannot be null");

        this.inputReader = inputReader;
        this.inputTimer = inputTimer;
    }

    public synchronized boolean hasNextNow() {
        return inputReader.hasNext();
    }

    public synchronized InputValue next() {
        inputTimer.restart();

        final InputValue v;
        if ( inputReader.hasNext() ) {
            v = new InputValue(inputReader.next());
        } else {
            v = new InputValue();
        }

        inputTimer.stop();

        return v;
    }

    /**
     * Helper class that contains a read value suitable for EOF marking in a BlockingQueue
     */
    class InputValue extends BlockingQueueValue<InputType> {
        private InputValue(InputType datum) { super(datum); }
        private InputValue() { }
    }
}
