package org.broadinstitute.sting.gatk.walkers.genotyper;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

public abstract class GenotypeLikelihoods {
    // precalculate these for performance (pow/log10 is expensive!)

    /**
     * SOLID data uses Q0 bases to represent reference-fixed bases -- they shouldn't be counted
     * during GL calculations.  If this field is true, Q0 bases will be removed in add().
     */
    protected boolean filterQ0Bases = true;

    protected static final double[] oneMinusData = new double[Byte.MAX_VALUE];
    protected static final double[] oneHalfMinusDataArachne = new double[Byte.MAX_VALUE];
    protected static final double[] oneHalfMinusData3Base = new double[Byte.MAX_VALUE];
    protected static final double log10Of1_3 = log10(1.0 / 3);
    protected static final double log10Of2_3 = log10(2.0 / 3);

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10, (qual / -10.0)));
        }
    }

    static {
        for (int qual = 0; qual < Byte.MAX_VALUE; qual++) {
            double e = pow(10, (qual / -10.0));
            oneHalfMinusDataArachne[qual] = log10(0.5 - e / 2.0);
            oneHalfMinusData3Base[qual] = log10(0.5 - e / 2.0 + e / 6.0);
            //System.out.printf("qual=%d, e=%f, oneHalfMinusDataArachne=%f, oneHalfMinusData3Base=%f%n", qual, e, oneHalfMinusDataArachne[qual], oneHalfMinusData3Base[qual]);
        }
    }

    public final static String[] genotypes = new String[10];
    static {
        genotypes[0] = "AA";
        genotypes[1] = "AC";
        genotypes[2] = "AG";
        genotypes[3] = "AT";
        genotypes[4] = "CC";
        genotypes[5] = "CG";
        genotypes[6] = "CT";
        genotypes[7] = "GG";
        genotypes[8] = "GT";
        genotypes[9] = "TT";
    }

    public static final double HUMAN_HETEROZYGOSITY = 1e-3;

    public static double[] computePriors(double h) {
        double[] pdbls = new double[3];
        pdbls[0] = 1.0 - (3.0 * h / 2.0);
        pdbls[1] = h;
        pdbls[2] = h / 2.0;
        return pdbls;
    }

    public double[] likelihoods;

    public abstract int add(char ref, char read, byte qual);
    public abstract void applyPrior(char ref);
}
