package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;

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
        }
    }

    public static final double HUMAN_HETEROZYGOSITY = 1e-3;

    /**
     * Returns homozygous-reference, heterozygous, and homozygous-non-ref probabilities given a heterozygosity
     * value, as elements 0, 1, and 2 of a double[], respectively
     *
     * @param h the heterozygosity [probability of a base being heterozygous]
     */
    public static double[] heterozygosity2DiploidProbabilities(double h) {
        if (MathUtils.isNegativeOrZero(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        double[] pdbls = new double[3];
        pdbls[0] = 1.0 - (3.0 * h / 2.0);
        pdbls[1] = h;
        pdbls[2] = h / 2.0;
        return pdbls;
    }

    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * @param ref
     * @param priorHomRef
     * @param priorHet
     * @param priorHomVar
     */
    public static double[] getGenotypePriors(char ref, double priorHomRef, double priorHet, double priorHomVar) {
        if ((priorHomRef + priorHet + priorHomVar) != 1) {
            throw new RuntimeException(String.format("Prior probabilities don't sum to one => %f, %f, %f", priorHomRef, priorHet, priorHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int nRefBases = Utils.countOccurrences(ref, g.toString());
            double log10POfG = 0.0;

            switch ( nRefBases ) {
                case 2: // hom-ref
                    log10POfG = Math.log10(priorHomRef);
                    break;
                case 0: // hom-nonref
                    //log10POfG = Math.log10(priorHomVar / 3);
                    log10POfG = Math.log10(priorHomVar);
                    break;
                default:
                    //log10POfG = Math.log10(priorHet / 6);
                    log10POfG = Math.log10(priorHet);
                    break;
            }

            priors[g.ordinal()] = log10POfG;
        }

        return priors;
    }

    public abstract int add(char ref, char read, byte qual);
}
