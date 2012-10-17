package org.broadinstitute.sting.utils.nanoScheduler;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.AutoFormattingTime;
import org.broadinstitute.sting.utils.SimpleTimer;

/**
 * Holds runtime profile (input, read, map) times as tracked by NanoScheduler
 *
 * User: depristo
 * Date: 9/10/12
 * Time: 8:31 PM
 */
public class NSRuntimeProfile {
    final SimpleTimer outsideSchedulerTimer = new SimpleTimer("outside");
    final SimpleTimer inputTimer = new SimpleTimer("input");
    final SimpleTimer mapTimer = new SimpleTimer("map");
    final SimpleTimer reduceTimer = new SimpleTimer("reduce");

    /**
     * Combine the elapsed time information from other with this profile
     *
     * @param other a non-null profile
     */
    public void combine(final NSRuntimeProfile other) {
        outsideSchedulerTimer.addElapsed(other.outsideSchedulerTimer);
        inputTimer.addElapsed(other.inputTimer);
        mapTimer.addElapsed(other.mapTimer);
        reduceTimer.addElapsed(other.reduceTimer);
    }

    /**
     * Print the runtime profiling to logger
     *
     * @param logger
     */
    public void log(final Logger logger) {
        log1(logger, "Input   time", inputTimer);
        log1(logger, "Map     time", mapTimer);
        log1(logger, "Reduce  time", reduceTimer);
        log1(logger, "Outside time", outsideSchedulerTimer);
    }

    /**
     * @return the total runtime for all functions of this nano scheduler
     */
    //@Ensures("result >= 0.0")
    public double totalRuntimeInSeconds() {
        return inputTimer.getElapsedTime()
                + mapTimer.getElapsedTime()
                + reduceTimer.getElapsedTime()
                + outsideSchedulerTimer.getElapsedTime();
    }

    /**
     * Print to logger.info timing information from timer, with name label
     *
     * @param label the name of the timer to display.  Should be human readable
     * @param timer the timer whose elapsed time we will display
     */
    //@Requires({"label != null", "timer != null"})
    private void log1(final Logger logger, final String label, final SimpleTimer timer) {
        final double myTimeInSec = timer.getElapsedTime();
        final double myTimePercent = myTimeInSec / totalRuntimeInSeconds() * 100;
        logger.info(String.format("%s: %s (%5.2f%%)", label, new AutoFormattingTime(myTimeInSec), myTimePercent));
    }
}
