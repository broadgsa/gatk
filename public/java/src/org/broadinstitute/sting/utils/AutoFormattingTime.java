package org.broadinstitute.sting.utils;

/**
 * Simple utility class that makes it convenient to print unit adjusted times
 */
public class AutoFormattingTime {
    private final int width; // for format
    private final int precision;      // for format

    double timeInSeconds;           // in Seconds
    private final String formatString;

    public AutoFormattingTime(double timeInSeconds, final int width, int precision) {
        this.width = width;
        this.timeInSeconds = timeInSeconds;
        this.precision = precision;
        this.formatString = "%" + width + "." + precision + "f %s";
    }

    public AutoFormattingTime(double timeInSeconds, int precision) {
        this(timeInSeconds, 6, precision);
    }

    public AutoFormattingTime(double timeInSeconds) {
        this(timeInSeconds, 1);
    }

    public double getTimeInSeconds() {
        return timeInSeconds;
    }

    /**
     * @return the precision (a la format's %WIDTH.PERCISIONf)
     */
    public int getWidth() {
        return width;
    }

    /**
     * @return the precision (a la format's %WIDTH.PERCISIONf)
     */
    public int getPrecision() {
        return precision;
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

        return String.format(formatString, unitTime, unit);
    }
}
