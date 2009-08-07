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
    protected static final double log10Of1_3 = log10(1.0 / 3.0);
    
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

    public abstract int add(char ref, char read, byte qual);
}
