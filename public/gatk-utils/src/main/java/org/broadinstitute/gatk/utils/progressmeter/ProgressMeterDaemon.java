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

package org.broadinstitute.gatk.utils.progressmeter;

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
            meter.updateElapsedTimeInNanoseconds();
            try {
                Thread.sleep(getPollFrequencyMilliseconds());
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
