package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMUtils;

/**
 * QualityUtils is a static class (no instantiation allowed!) with some utility methods for manipulating
 * quality scores.
 *
 * @author Kiran Garimella
 */
public class QualityUtils {
    public final static byte MAX_RECALIBRATED_Q_SCORE = SAMUtils.MAX_PHRED_SCORE;
    public final static byte MAX_QUAL_SCORE = SAMUtils.MAX_PHRED_SCORE;
    public final static double ERROR_RATE_OF_MAX_QUAL_SCORE = qualToErrorProbRaw(MAX_QUAL_SCORE);

    public final static double MIN_REASONABLE_ERROR = 0.0001;
    public final static byte MAX_REASONABLE_Q_SCORE = 60;  // quals above this value are extremely suspicious
    public final static byte MIN_USABLE_Q_SCORE = 6;
    public final static int MAPPING_QUALITY_UNAVAILABLE = 255;

    private static double qualToErrorProbCache[] = new double[256];
    static {
        for (int i = 0; i < 256; i++) qualToErrorProbCache[i] = qualToErrorProbRaw(i);
    }

    private static double qualToErrorProbLog10Cache[] = new double[256];
    static {
        for (int i = 0; i < 256; i++) qualToErrorProbLog10Cache[i] = qualToErrorProbLog10Raw(i);
    }

    private static double qualToProbLog10Cache[] = new double[256];
    static {
        for (int i = 0; i < 256; i++) qualToProbLog10Cache[i] = qualToProbLog10Raw(i);
    }

    /**
     * Private constructor.  No instantiating this class!
     */
    private QualityUtils() {}

    /**
     * Convert a quality score to a probability.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0-255)
     * @return a probability (0.0-1.0)
     */
    static public double qualToProb(byte qual) {
        return 1.0 - qualToErrorProb(qual);
    }

    static public double qualToProb(double qual) {
        return 1.0 - Math.pow(10.0, qual/(-10.0));
    }

    static private double qualToProbLog10Raw(int qual) {
        return Math.log10(1.0 - qualToErrorProbRaw(qual));
    }

    static public double qualToProbLog10(byte qual) {
        return qualToProbLog10Cache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    /**
     * Convert a quality score to a probability of error.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0 - 255)
     * @return a probability (0.0 - 1.0)
     */
    static private double qualToErrorProbRaw(int qual) {
        return qualToErrorProb((double) qual);
    }

    public static double qualToErrorProb(final double qual) {
        return Math.pow(10.0, ((double) qual)/-10.0);
    }


    static public double qualToErrorProb(byte qual) {
        return qualToErrorProbCache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    static private double qualToErrorProbLog10Raw(int qual) {
        return ((double) qual)/-10.0;
    }

    static public double qualToErrorProbLog10(byte qual) {
        return qualToErrorProbLog10Cache[(int)qual & 0xff]; // Map: 127 -> 127; -128 -> 128; -1 -> 255; etc.
    }

    /**
     * Convert a probability to a quality score.  Note, this is capped at Q40.
     *
     * @param prob a probability (0.0-1.0)
     * @return a quality score (0-40)
     */
    static public byte probToQual(double prob) {
        return probToQual(prob, MIN_REASONABLE_ERROR);
        //return (byte) Math.round(-10.0*Math.log10(1.0 - prob + 0.0001));
    }

    /**
     * Convert a probability to a quality score.  Note, this is capped at a quality score which is determined by _eps_.
     *
     * @param prob a probability (0.0-1.0)
     * @param eps min probabilty allowed (0.0-1.0)
     * @return a quality score (0-255)
     */
    static public byte probToQual(double prob, double eps) {
        double lp = Math.round(-10.0*Math.log10(1.0 - prob + eps));
        //System.out.printf("LP is %f, byte is %d%n", lp, b);
        return boundQual((int)lp);
    }

    static public double phredScaleCorrectRate(double trueRate) {
        return phredScaleErrorRate(1-trueRate);
    }

    static public double phredScaleErrorRate(double errorRate) {
        return Math.abs(-10.0*Math.log10(errorRate));
    }

    /**
     * Return a quality score, capped at max qual.
     *
     * @param qual  the uncapped quality score
     * @return the capped quality score
     */
    static public byte boundQual(int qual) {
        return boundQual(qual, MAX_QUAL_SCORE);
    }

    /**
     * Returns an integer quality score bounded by 1 - maxQual.
     *
     * @param qual    the quality score
     * @param maxQual the maximum quality
     * @return the integer betwen 1 and maxqual.
     */
    static public byte boundQual(int qual, byte maxQual) {
        return (byte) Math.max(Math.min(qual, maxQual), 1);
    }
}
