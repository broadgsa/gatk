package org.broadinstitute.sting.utils;

import java.io.PrintStream;
import com.google.java.contract.*;

/**
 * A useful simple system for timing code.  This code is not thread safe!
 *
 * User: depristo
 * Date: Dec 10, 2010
 * Time: 9:07:44 AM
 */
@Invariant({
        "elapsed >= 0",
        "startTime >= 0",
        "name != null",
        "! running || startTime > 0"})
public class SimpleTimer {
    final private String name;
    private long elapsed = 0l;
    private long startTime = 0l;
    boolean running = false;

    /**
     * Creates an anonymous simple timer
     */
    @Ensures("name != null && name.equals(\"Anonymous\")")
    public SimpleTimer() {
        this("Anonymous");
    }

    /**
     * Creates a simple timer named name
     * @param name of the timer, must not be null
     */
    @Requires("name != null")
    @Ensures("this.name != null && this.name.equals(name)")
    public SimpleTimer(String name) {
        this.name = name;
    }

    /**
     * @return the name associated with this timer
     */
    @Ensures("result != null")
    public String getName() {
        return name;
    }

    /**
     * Starts the timer running, and sets the elapsed time to 0.  This is equivalent to
     * resetting the time to have no history at all.
     *
     * @return this object, for programming convenience
     */
    @Requires("running == false")
    @Ensures({"result != null", "elapsed == 0l"})
    public SimpleTimer start() {
        elapsed = 0l;
        restart();
        return this;
    }

    /**
     * Starts the timer running, without reseting the elapsed time.  This function may be
     * called without first calling start().  The only difference between start and restart
     * is that start resets the elapsed time, while restart does not.
     *
     * @return this object, for programming convenience
     */
    @Requires("running == false")
    @Ensures("result != null")
    public SimpleTimer restart() {
        running = true;
        startTime = currentTime();
        return this;
    }

    /**
     * @return is this timer running?
     */
    public boolean isRunning() {
        return running;
    }

    /**
     * @return A convenience function to obtain the current time in milliseconds from this timer
     */
    public long currentTime() {
        return System.currentTimeMillis();
    }

    /**
     * Stops the timer.  Increases the elapsed time by difference between start and now.  The
     * timer must be running in order to call stop
     *
     * @return this object, for programming convenience
     */
    @Requires("running == true")
    @Ensures({"result != null", "elapsed >= old(elapsed)", "running == false"})
    public SimpleTimer stop() {
        running = false;
        elapsed += currentTime() - startTime;
        return this;
    }

    /**
     * Returns the total elapsed time of all start/stops of this timer.  If the timer is currently
     * running, includes the difference from currentTime() and the start as well
     *
     * @return this time, in seconds
     */
    @Ensures({
            "result >= (elapsed/1000.0)",
            "result >= 0"})
    public double getElapsedTime() {
        return (running ? (currentTime() - startTime + elapsed) : elapsed) / 1000.0;
    }


    public void printElapsedTime(PrintStream out) {
        out.printf("SimpleTimer %s: %.2f%n", name, getElapsedTime());
    }
}
