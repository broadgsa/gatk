package org.broadinstitute.sting.utils;

/**
 * Simple utility class that makes it convenient to print unit adjusted times
 */
public class AutoFormattingTime {
    double timeInSeconds;           // in Seconds
    int precision;      // for format

    public AutoFormattingTime(double timeInSeconds, int precision) {
        this.timeInSeconds = timeInSeconds;
        this.precision = precision;
    }

    public AutoFormattingTime(double timeInSeconds) {
        this(timeInSeconds, 1);
    }

    public double getTimeInSeconds() {
        return timeInSeconds;
    }

    /**
     * Instead of 10000 s, returns 2.8 hours
     * @return
     */
    public String toString() {
        double unitTime = timeInSeconds;
        String unit = "s";

        if ( timeInSeconds > 120 ) {
            unitTime = timeInSeconds / 60; // minutes
            unit = "m";

            if ( unitTime > 120 ) {
                unitTime /= 60; // hours
                unit = "h";

                if ( unitTime > 100 ) {
                    unitTime /= 24; // days
                    unit = "d";

                    if ( unitTime > 20 ) {
                        unitTime /= 7; // days
                        unit = "w";
                    }
                }
            }
        }

        return String.format("%6."+precision+"f %s", unitTime, unit);
    }
}
