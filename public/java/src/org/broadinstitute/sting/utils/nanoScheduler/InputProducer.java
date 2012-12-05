package org.broadinstitute.sting.utils.nanoScheduler;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MultiThreadedErrorTracker;

import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CountDownLatch;

/**
 * Producer Thread that reads input values from an inputReads and puts them into an output queue
 */
class InputProducer<InputType> implements Runnable {
    private final static Logger logger = Logger.getLogger(InputProducer.class);

    /**
     * The iterator we are using to get data from
     */
    final Iterator<InputType> inputReader;

    /**
     * Where we put our input values for consumption
     */
    final BlockingQueue<InputValue> outputQueue;

    final MultiThreadedErrorTracker errorTracker;

    /**
     * Have we read the last value from inputReader?
     *
     * Must be a local variable, as inputReader.hasNext() can actually end up doing a lot
     * of work, and the method getNumInputValues() is supposed to be called not in the
     * thread executing the reading of values but in the thread enqueuing results
     */
    boolean readLastValue = false;

    int nRead = 0;
    int inputID = -1;

    /**
     * A latch used to block threads that want to start up only when all of the values
     * in inputReader have been read by the thread executing run()
     */
    final CountDownLatch latch = new CountDownLatch(1);

    public InputProducer(final Iterator<InputType> inputReader,
                         final MultiThreadedErrorTracker errorTracker,
                         final BlockingQueue<InputValue> outputQueue) {
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        if ( errorTracker == null ) throw new IllegalArgumentException("errorTracker cannot be null");
        if ( outputQueue == null ) throw new IllegalArgumentException("OutputQueue cannot be null");

        this.inputReader = inputReader;
        this.errorTracker = errorTracker;
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
        if ( ! inputReader.hasNext() ) {
            // we are done, mark ourselves as such and return null
            readLastValue = true;
            return null;
        } else {
            // get the next value, and return it
            final InputType input = inputReader.next();
            if ( input == null )
                throw new IllegalStateException("inputReader.next() returned a null value, breaking our contract");
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
                final InputType value = readNextItem();

                if ( value == null ) {
                    if ( ! readLastValue )
                        throw new IllegalStateException("value == null but readLastValue is false!");

                    // add the EOF object so our consumer knows we are done in all inputs
                    // note that we do not increase inputID here, so that variable indicates the ID
                    // of the last real value read from the queue
                    outputQueue.put(new InputValue(inputID + 1));
                    break;
                } else {
                    // add the actual value to the outputQueue
                    outputQueue.put(new InputValue(++inputID, value));
                }
            }

            latch.countDown();
        } catch (Throwable ex) {
            errorTracker.notifyOfError(ex);
        } finally {
//            logger.info("Exiting input thread readLastValue = " + readLastValue);
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
     *
     * This class also contains an ID, an integer incrementing from 0 to N, for N total
     * values in the input stream.  This ID indicates which element in the element stream this
     * InputValue corresponds to.  Necessary for tracking and ordering results by input position.
     *
     * Note that EOF markers have IDs > N, and ID values >> N can occur if many EOF markers
     * are enqueued in the outputQueue.
     */
    class InputValue extends EOFMarkedValue<InputType> {
        final int id;

        private InputValue(final int id, InputType datum) {
            super(datum);
            if ( id < 0 ) throw new IllegalArgumentException("id must be >= 0");
            this.id = id;
        }
        private InputValue(final int id) {
            super();
            if ( id < 0 ) throw new IllegalArgumentException("id must be >= 0");
            this.id = id;
        }

        /**
         * Returns the ID of this input marker
         * @return id >= 0
         */
        public int getId() {
            return id;
        }

        /**
         * Create another EOF marker with ID + 1 to this one.
         *
         * Useful in the case where we need to enqueue another EOF marker for future jobs and we
         * want them to have a meaningful ID, one greater than the last one.
         *
         * @return ID
         */
        //@Ensures({"result.isEOFMarker()", "result.getId() == getId() + 1"})
        public InputValue nextEOF() {
            if ( ! isEOFMarker() )
                throw new IllegalArgumentException("Cannot request next EOF marker for non-EOF marker InputValue");
            return new InputValue(getId() + 1);
        }
    }
}
