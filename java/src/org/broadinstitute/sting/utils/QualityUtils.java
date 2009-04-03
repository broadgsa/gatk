package org.broadinstitute.sting.utils;

public class QualityUtils {
    private QualityUtils() {}
    
    static public double qualToProb(byte qual) {
        return 1.0 - Math.pow(10.0, ((double) qual)/-10.0);
    }

    static public byte probToQual(double prob) {
        return (byte) Math.round(-10.0*Math.log10(1.0 - prob + 0.0001));
    }

    static public byte baseAndProbToCompressedQuality(int baseIndex, double prob) {
        byte compressedQual = (byte) baseIndex;
        byte cprob = (byte) (100.0*prob);
        byte qualmask = (byte) 252;
        compressedQual += ((cprob << 2) & qualmask);

        return compressedQual;
    }
}
