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

package org.broadinstitute.gatk.utils;

import java.util.concurrent.TimeUnit;

/**
 * Conveniently print a time with an automatically determined time unit
 *
 * For example, if the amount of time is 10^6 seconds, instead of printing
 * out 10^6 seconds, prints out 11.57 days instead.
 *
 * Dynamically uses time units:
 *
 *   - seconds: s
 *   - minutes: m
 *   - hours  : h
 *   - days   : d
 *   - weeks  : w
 *
 * @author depristo
 * @since 2009
 */
public class AutoFormattingTime {
    private static final double NANOSECONDS_PER_SECOND = 1e9;

    /**
     * Width a la format's %WIDTH.PERCISIONf
     */
    private final int width; // for format

    /**
     * Precision a la format's %WIDTH.PERCISIONf
     */
    private final int precision;      // for format

    /**
     * The elapsed time in nanoseconds
     */
    private final long nanoTime;

    /**
     * Create a new autoformatting time with elapsed time nanoTime in nanoseconds
     * @param nanoTime the elapsed time in nanoseconds
     * @param width the width >= 0 (a la format's %WIDTH.PERCISIONf) to use to display the format, or -1 if none is required
     * @param precision the precision to display the time at.  Must be >= 0;
     */
    public AutoFormattingTime(final long nanoTime, final int width, int precision) {
        if ( width < -1 ) throw new IllegalArgumentException("Width " + width + " must be >= -1");
        if ( precision < 0 ) throw new IllegalArgumentException("Precision " + precision + " must be >= 0");

        this.width = width;
        this.nanoTime = nanoTime;
        this.precision = precision;
    }

    /**
     * @see #AutoFormattingTime(long, int, int) but with default width and precision
     * @param nanoTime
     */
    public AutoFormattingTime(final long nanoTime) {
        this(nanoTime, 6, 1);
    }

    /**
     * @see #AutoFormattingTime(long, int, int) but with time specificied as a double in seconds
     */
    public AutoFormattingTime(final double timeInSeconds, final int width, final int precision) {
        this(secondsToNano(timeInSeconds), width, precision);
    }

    /**
     * @see #AutoFormattingTime(long) but with time specificied as a double in seconds
     */
    public AutoFormattingTime(double timeInSeconds) {
        this(timeInSeconds, 6, 1);
    }

    /**
     * Precomputed format string suitable for string.format with the required width and precision
     */
    private String getFormatString() {
        final StringBuilder b = new StringBuilder("%");
        if ( width != -1 )
            b.append(width);
        b.append(".").append(precision).append("f %s");
        return b.toString();
    }

    /**
     * Get the time associated with this object in nanoseconds
     * @return the time in nanoseconds
     */
    public long getTimeInNanoSeconds() {
        return nanoTime;
    }

    /**
     * Get the time associated with this object in seconds, as a double
     * @return time in seconds as a double
     */
    public double getTimeInSeconds() {
        return TimeUnit.NANOSECONDS.toSeconds(getTimeInNanoSeconds());
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
     * Get a string representation of this time, automatically converting the time
     * to a human readable unit with width and precision provided during construction
     * @return a non-null string
     */
    public String toString() {
        double unitTime = getTimeInSeconds();
        String unit = "s";

        if ( unitTime > 120 ) {
            unitTime /= 60; // minutes
            unit = "m";

            if ( unitTime > 120 ) {
                unitTime /= 60; // hours
                unit = "h";

                if ( unitTime > 100 ) {
                    unitTime /= 24; // days
                    unit = "d";

                    if ( unitTime > 20 ) {
                        unitTime /= 7; // weeks
                        unit = "w";
                    }
                }
            }
        }

        return String.format(getFormatString(), unitTime, unit);
    }


    /**
     * Convert a time in seconds as a double into nanoseconds as a long
     * @param timeInSeconds an elapsed time in seconds, as a double
     * @return an equivalent value in nanoseconds as a long
     */
    private static long secondsToNano(final double timeInSeconds) {
        return (long)(NANOSECONDS_PER_SECOND * timeInSeconds);
    }

}
