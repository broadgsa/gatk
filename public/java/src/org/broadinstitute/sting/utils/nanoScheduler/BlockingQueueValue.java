package org.broadinstitute.sting.utils.nanoScheduler;

import com.google.java.contract.Invariant;

/**
 * Wrapper to hold data for a blocking queue, distinguishing an EOF marker from a real object
 *
 * The only way to tell in a consumer thread that a blocking queue has no more data ever
 * coming down the pipe is to pass in a "poison" or EOF object.  This class provides
 * a generic capacity for that...
 *
 * The use case looks like this:
 *
 * BlockingQueue q
 * producer:
 *   while ( x has items )
 *      q.put(new BlockingQueueValue(x))
 *   q.put(new BlockingQueueValue())
 *
 * Consumer:
 *   while ( true )
 *       value = q.take()
 *       if ( value.isLast() )
 *          break
 *       else
 *          do something useful with value
 *
 *
 * User: depristo
 * Date: 9/6/12
 * Time: 3:08 PM
 */
@Invariant("! isLast || value == null")
class BlockingQueueValue<T> {
    /**
     * True if this is the EOF marker object
     */
    final private boolean isLast;

    /**
     * Our value, if we aren't the EOF marker
     */
    final private T value;

    /**
     * Create a new BlockingQueueValue containing a real value, where last is false
     * @param value
     */
    BlockingQueueValue(final T value) {
        isLast = false;
        this.value = value;
    }

    /**
     * Create a new BlockingQueueValue that is the last item
     */
    BlockingQueueValue() {
        isLast = true;
        this.value = null;
    }

    /**
     * Is this the EOF marker?
     *
     * @return true if so, else false
     */
    public boolean isLast() {
        return isLast;
    }

    /**
     * Get the value held by this BlockingQueueValue
     *
     * @return the value
     * @throws IllegalStateException if this is the last item
     */
    public T getValue() {
        if ( isLast() )
            throw new IllegalStateException("Cannot get value for last object");
        return value;
    }
}
