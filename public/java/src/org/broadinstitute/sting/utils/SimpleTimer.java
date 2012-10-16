package org.broadinstitute.sting.utils;


import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.concurrent.TimeUnit;

/**
 * A useful simple system for timing code with nano second resolution
 *
 * Note that this code is not thread-safe.  If you have a single timer
 * being started and stopped by multiple threads you will need to protect the
 * calls to avoid meaningless results of having multiple starts and stops
 * called sequentially.
 *
 * User: depristo
 * Date: Dec 10, 2010
 * Time: 9:07:44 AM
 */
public class SimpleTimer {
    protected static final double NANO_TO_SECOND_DOUBLE = 1.0 / TimeUnit.SECONDS.toNanos(1);
    private final String name;

    /**
     * The elapsedTimeNano time in nanoSeconds of this timer.  The elapsedTimeNano time is the
     * sum of times between starts/restrats and stops.
     */
    private long elapsedTimeNano = 0l;

    /**
     * The start time of the last start/restart in nanoSeconds
     */
    private long startTimeNano = 0l;

    /**
     * Is this timer currently running (i.e., the last call was start/restart)
     */
    private boolean running = false;

    /**
     * Creates an anonymous simple timer
     */
    public SimpleTimer() {
        this("Anonymous");
    }

    /**
     * Creates a simple timer named name
     * @param name of the timer, must not be null
     */
    public SimpleTimer(final String name) {
        if ( name == null ) throw new IllegalArgumentException("SimpleTimer name cannot be null");
        this.name = name;
    }

    /**
     * @return the name associated with this timer
     */
    public synchronized String getName() {
        return name;
    }

    /**
     * Starts the timer running, and sets the elapsedTimeNano time to 0.  This is equivalent to
     * resetting the time to have no history at all.
     *
     * @return this object, for programming convenience
     */
    @Ensures("elapsedTimeNano == 0l")
    public synchronized SimpleTimer start() {
        elapsedTimeNano = 0l;
        return restart();
    }

    /**
     * Starts the timer running, without resetting the elapsedTimeNano time.  This function may be
     * called without first calling start().  The only difference between start and restart
     * is that start resets the elapsedTimeNano time, while restart does not.
     *
     * @return this object, for programming convenience
     */
    public synchronized SimpleTimer restart() {
        running = true;
        startTimeNano = currentTimeNano();
        return this;
    }

    /**
     * @return is this timer running?
     */
    public synchronized boolean isRunning() {
        return running;
    }

    /**
     * @return A convenience function to obtain the current time in milliseconds from this timer
     */
    public long currentTime() {
        return System.currentTimeMillis();
    }

    /**
     * @return A convenience function to obtain the current time in nanoSeconds from this timer
     */
    public long currentTimeNano() {
        return System.nanoTime();
    }

    /**
     * Stops the timer.  Increases the elapsedTimeNano time by difference between start and now.
     *
     * It's ok to call stop on a timer that's not running.  It has no effect on the timer.
     *
     * @return this object, for programming convenience
     */
    @Requires("startTimeNano != 0l")
    public synchronized SimpleTimer stop() {
        if ( running ) {
            running = false;
            elapsedTimeNano += currentTimeNano() - startTimeNano;
        }
        return this;
    }

    /**
     * Returns the total elapsedTimeNano time of all start/stops of this timer.  If the timer is currently
     * running, includes the difference from currentTime() and the start as well
     *
     * @return this time, in seconds
     */
    public synchronized double getElapsedTime() {
        return nanoToSecondsAsDouble(getElapsedTimeNano());
    }

    protected static double nanoToSecondsAsDouble(final long nano) {
        return nano * NANO_TO_SECOND_DOUBLE;
    }

    /**
     * @see #getElapsedTime() but returns the result in nanoseconds
     *
     * @return the elapsed time in nanoseconds
     */
    public synchronized long getElapsedTimeNano() {
        return running ? (currentTimeNano() - startTimeNano + elapsedTimeNano) : elapsedTimeNano;
    }

    /**
     * Add the elapsed time from toAdd to this elapsed time
     *
     * @param toAdd the timer whose elapsed time we want to add to this timer
     */
    public synchronized void addElapsed(final SimpleTimer toAdd) {
        elapsedTimeNano += toAdd.getElapsedTimeNano();
    }
}
