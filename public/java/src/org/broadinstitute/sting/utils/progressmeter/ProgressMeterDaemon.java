package org.broadinstitute.sting.utils.progressmeter;

/**
 * Daemon thread that periodically prints the progress of the progress meter
 *
 * User: depristo
 * Date: 12/4/12
 * Time: 9:16 PM
 */
public final class ProgressMeterDaemon extends Thread {
    public final static long DEFAULT_POLL_FREQUENCY_MILLISECONDS = 10 * 1000;

    /**
     * How frequently should we poll and print progress?
     */
    private final long pollFrequencyMilliseconds;

    /**
     * How long are we waiting between print progress calls are issued?
     * @return the time in milliseconds between progress meter calls
     */
    private long getPollFrequencyMilliseconds() {
        return pollFrequencyMilliseconds;
    }

    /**
     * Are we to continue periodically printing status, or should we shut down?
     */
    boolean done = false;

    /**
     * The meter we will call print on
     */
    final ProgressMeter meter;

    /**
     * Create a new ProgressMeterDaemon printing progress for meter
     * @param meter the progress meter to print progress of
     */
    public ProgressMeterDaemon(final ProgressMeter meter, final long pollFrequencyMilliseconds) {
        if ( meter == null ) throw new IllegalArgumentException("meter cannot be null");
        if ( pollFrequencyMilliseconds <= 0 ) throw new IllegalArgumentException("pollFrequencyMilliseconds must be greater than 0 but got " + pollFrequencyMilliseconds);

        this.meter = meter;
        this.pollFrequencyMilliseconds = pollFrequencyMilliseconds;
        setDaemon(true);
        setName("ProgressMeterDaemon");
    }

    public ProgressMeterDaemon(final ProgressMeter meter) {
        this(meter, DEFAULT_POLL_FREQUENCY_MILLISECONDS);
    }

    /**
     * Tells this daemon thread to shutdown at the next opportunity, as the progress
     * metering is complete.
     */
    public final void done() {
        this.done = true;
    }

    /**
     * Is this daemon thread done?
     * @return true if done, false otherwise
     */
    public boolean isDone() {
        return done;
    }

    /**
     * Start up the ProgressMeterDaemon, polling every tens of seconds to print, if
     * necessary, the provided progress meter.  Never exits until the JVM is complete,
     * or done() is called, as the thread is a daemon thread
     */
    public void run() {
        while (! done) {
            meter.printProgress(false);
            try {
                Thread.sleep(getPollFrequencyMilliseconds());
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
