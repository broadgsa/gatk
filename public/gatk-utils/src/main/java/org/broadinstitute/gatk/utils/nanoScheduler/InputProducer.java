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

import org.apache.log4j.Logger;

import java.util.Iterator;
import java.util.concurrent.CountDownLatch;

/**
 * Helper class that allows multiple threads to reads input values from
 * an iterator, and track the number of items read from that iterator.
 */
class InputProducer<InputType> {
    private final static Logger logger = Logger.getLogger(InputProducer.class);

    /**
     * The iterator we are using to get data from
     */
    final Iterator<InputType> inputReader;

    /**
     * Have we read the last value from inputReader?
     *
     * Must be a local variable, as inputReader.hasNext() can actually end up doing a lot
     * of work, and the method getNumInputValues() is supposed to be called not in the
     * thread executing the reading of values but in the thread enqueuing results
     */
    boolean readLastValue = false;

    /**
     * Once we've readLastValue, lastValue contains a continually
     * updating InputValue where EOF is true.  It's not necessarily
     * a single value, as each read updates lastValue with the
     * next EOF marker
     */
    private InputValue lastValue = null;

    int nRead = 0;
    int inputID = -1;

    public InputProducer(final Iterator<InputType> inputReader) {
        if ( inputReader == null ) throw new IllegalArgumentException("inputReader cannot be null");
        this.inputReader = inputReader;
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
     */
    private synchronized InputType readNextItem() {
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
     * Are there currently more values in the iterator?
     *
     * Note the word currently.  It's possible that some already submitted
     * job will read a value from this InputProvider, so in some sense
     * there are no more values and in the future there'll be no next
     * value.  That said, once this returns false it means that all
     * of the possible values have been read
     *
     * @return true if a future call to next might return a non-EOF value, false if
     *         the underlying iterator is definitely empty
     */
    public synchronized boolean hasNext() {
        return ! allInputsHaveBeenRead();
    }

    /**
     * Get the next InputValue from this producer.  The next value is
     * either (1) the next value from the iterator, in which case the
     * the return value is an InputValue containing that value, or (2)
     * an InputValue with the EOF marker, indicating that the underlying
     * iterator has been exhausted.
     *
     * This function never fails -- it can be called endlessly and
     * while the underlying iterator has values it returns them, and then
     * it returns a succession of EOF marking input values.
     *
     * @return an InputValue containing the next value in the underlying
     *         iterator, or one with EOF marker, if the iterator is exhausted
     */
    public synchronized InputValue next() {
        if ( readLastValue ) {
            // we read the last value, so our value is the next
            // EOF marker based on the last value.  Make sure to
            // update the last value so the markers keep incrementing
            // their job ids
            lastValue = lastValue.nextEOF();
            return lastValue;
        } else {
            final InputType value = readNextItem();

            if ( value == null ) {
                if ( ! readLastValue )
                    throw new IllegalStateException("value == null but readLastValue is false!");

                // add the EOF object so our consumer knows we are done in all inputs
                // note that we do not increase inputID here, so that variable indicates the ID
                // of the last real value read from the queue
                lastValue = new InputValue(inputID + 1);
                return lastValue;
            } else {
                // add the actual value to the outputQueue
                return new InputValue(++inputID, value);
            }
        }
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
