package org.broadinstitute.sting.utils.progressmeter;

/**
 * Daemon thread that periodically prints the progress of the progress meter
 *
 * User: depristo
 * Date: 12/4/12
 * Time: 9:16 PM
 */
public final class ProgressMeterDaemon extends Thread {
    /**
     * How frequently should we poll and print progress?
     */
    private final static long POLL_FREQUENCY_MILLISECONDS = 10 * 1000;

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
    public ProgressMeterDaemon(final ProgressMeter meter) {
        if ( meter == null ) throw new IllegalArgumentException("meter cannot be null");
        this.meter = meter;
        setDaemon(true);
        setName("ProgressMeterDaemon");
    }

    /**
     * Tells this daemon thread to shutdown at the next opportunity, as the progress
     * metering is complete.
     */
    public final void done() {
        this.done = true;
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
                Thread.sleep(POLL_FREQUENCY_MILLISECONDS);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
