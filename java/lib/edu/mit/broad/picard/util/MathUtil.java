package edu.mit.broad.picard.util;

/**
 * General math utilities
 *
 * @author Tim Fennell
 */
public class MathUtil {
    /** Calculated the mean of an array of doubles. */
    public static double mean(double[] in, int start, int stop) {
        double total = 0;
        for (int i=start; i<stop; ++i) {
            total += in[i];
        }

        return total / (stop-start);
    }

    /** Calculated the standard deviation of an array of doubles. */
    public static double stddev(double[] in, int start, int length) {
        return stddev(in, start, length, mean(in, start, length));
    }

    /** Calculated the standard deviation of an array of doubles. */
    public static double stddev(double[] in, int start, int stop, double mean) {
        double total = 0;
        for (int i=start; i<stop; ++i) {
            total += (in[i] * in[i]);
        }

        return Math.sqrt((total / (stop-start)) - (mean*mean));
    }
}
