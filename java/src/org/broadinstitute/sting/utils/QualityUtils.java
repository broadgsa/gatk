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
        return 1.0 - qualToErrorProb(qual);
    }

    static public double qualToProb(int qual) {
        return 1.0 - Math.pow(10.0, ((double) qual)/-10.0);
    }

    /**
     * Convert a quality score to a probability of error.  This is the Phred-style
     * conversion, *not* the Illumina-style conversion (though asymptotically, they're the same).
     *
     * @param qual a quality score (0-40)
     * @return a probability (0.0-1.0)
     */
    static public double qualToErrorProb(byte qual) {
        return Math.pow(10.0, ((double) qual)/-10.0);
    }

    /**
     * Convert a probability to a quality score.  Note, this is capped at Q40.
     *
     * @param prob a probability (0.0-1.0)
     * @return a quality score (0-40)
     */
    static public byte probToQual(double prob) {
        return probToQual(prob, 0.0001);
        //return (byte) Math.round(-10.0*Math.log10(1.0 - prob + 0.0001));
    }

    /**
     * Convert a probability to a quality score.  Note, this is capped at Q40.
     *
     * @param prob a probability (0.0-1.0)
     * @param eps min probabilty allowed (0.0-1.0)
     * @return a quality score (0-255)
     */
    static public byte probToQual(double prob, double eps) {
        double lp = Math.round(-10.0*Math.log10(1.0 - prob + eps));
        byte b = boundQual((int)lp);
        //System.out.printf("LP is %f, byte is %d%n", lp, b);
        return b;
    }

    /**
     * Return a quality score, capped at 63.
     *
     * @param qual  the uncapped quality score
     * @return the capped quality score
     */
    static public byte boundQual(int qual) {
        return (byte) Math.min(qual, 63);        
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
        byte compressedQual = 0;

        compressedQual = (byte) baseIndex;

        byte cprob = (byte) (100.0*prob);
        byte qualmask = (byte) 252;
        compressedQual += ((cprob << 2) & qualmask);
        
        return compressedQual;
    }

    /**
     * Compress a base and a log probabiliy difference (-10log10(p3/p2)) into
     * a single byte so that it can be output in a SAMRecord's SQ field.
     *
     * @param baseIndex  the base index
     * @param probdiff   the log probability difference between the secondary and tertiary bases (-10log10(p3/p2))
     * @return a byte containing the index and the log probability difference
     */
    static public byte baseAndProbDiffToCompressedQuality(int baseIndex, double probdiff) {
        byte compressedQual = 0;

        compressedQual = (byte) baseIndex;

        byte cprob = (byte) probdiff;
        byte qualmask = (byte) 252;
        compressedQual += ((cprob << 2) & qualmask);

        return compressedQual;
    }

    /**
     * From a compressed base, extract the base index (0:A, 1:C, 2:G, 3:T)
     *
     * @param compressedQual the compressed quality score, as returned by baseAndProbToCompressedQuality
     * @return base index
     */
    static public int compressedQualityToBaseIndex(byte compressedQual) {
        return (int) (compressedQual & 0x3);
    }

    /**
     * From a compressed base, extract the base probability
     *
     * @param compressedQual the compressed quality score, as returned by baseAndProbToCompressedQuality
     * @return the probability
     */
    static public double compressedQualityToProb(byte compressedQual) {
        // Because java natives are signed, extra care must be taken to avoid
        // shifting a 1 into the sign bit in the implicit promotion of 2 to an int.
        int x2 = ((int) compressedQual) & 0xff;
        x2 = (x2 >>> 2);

        return ((double) x2)/100.0;
    }

    /**
     * From a compressed base, extract the log probability difference between the secondary and tertiary bases.
     *
     * @param compressedQual the compressed quality score, as returned by baseAndProbDiffToCompressedQuality
     * @return the log probability difference (-10log10(p3/p2))
     */
    static public double compressedQualityToProbDiff(byte compressedQual) {
        // Because java natives are signed, extra care must be taken to avoid
        // shifting a 1 into the sign bit in the implicit promotion of 2 to an int.
        int x2 = ((int) compressedQual) & 0xff;
        x2 = (x2 >>> 2);

        return ((double) x2);
    }

    /**
     * Return the complement of a compressed quality
     *
     * @param compressedQual  the compressed quality score (as returned by baseAndProbToCompressedQuality)
     * @return the complementary compressed quality
     */
    static public byte complementCompressedQuality(byte compressedQual) {
        int baseIndex = compressedQualityToBaseIndex(compressedQual);
        double prob = compressedQualityToProb(compressedQual);

        return baseAndProbToCompressedQuality(BaseUtils.complementIndex(baseIndex), prob);
    }

    /**
     * Return the reverse complement of a byte array of compressed qualities
     *
     * @param compressedQuals  a byte array of compressed quality scores
     * @return the reverse complement of the byte array
     */
    static public byte[] reverseComplementCompressedQualityArray(byte[] compressedQuals) {
        byte[] rcCompressedQuals = new byte[compressedQuals.length];

        for (int pos = 0; pos < compressedQuals.length; pos++) {
            rcCompressedQuals[compressedQuals.length - pos - 1] = complementCompressedQuality(compressedQuals[pos]);
        }

        return rcCompressedQuals;
    }
}
