package org.broadinstitute.sting.utils;

/**
 * QualityUtils is a static class (no instantiation allowed!) with some utility methods for manipulating
 * quality scores.
 *
 * @author Kiran Garimella
 */
public class QualityUtils {
    /**
     * Private constructor.  No instantiating this class!
     */
    private QualityUtils() {}

    /**
     * Convert a quality score to a probability.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0-40)
     * @return a probability (0.0-1.0)
     */
    static public double qualToProb(byte qual) {
        return 1.0 - Math.pow(10.0, ((double) qual)/-10.0);
    }

    /**
     * Convert a probability to a quality score.  Note, this is capped at Q40.
     *
     * @param prob a probability (0.0-1.0)
     * @return a quality score (0-40)
     */
    static public byte probToQual(double prob) {
        return (byte) Math.round(-10.0*Math.log10(1.0 - prob + 0.0001));
    }

    /**
     * Compress a base and a probability into a single byte so that it can be output in a SAMRecord's SQ field.
     * Note: the highest probability this function can encode is 64%, so this function should only never be used on the best base hypothesis.
     * Another note: the probability encoded here gets rounded to the nearest 1%.
     *
     * @param baseIndex the base index
     * @param prob      the base probability
     * @return a byte containing the index and the probability
     */
    static public byte baseAndProbToCompressedQuality(int baseIndex, double prob) {
        byte compressedQual = (byte) baseIndex;
        byte cprob = (byte) (100.0*prob);
        byte qualmask = (byte) 252;
        compressedQual += ((cprob << 2) & qualmask);

        return compressedQual;
    }
}
