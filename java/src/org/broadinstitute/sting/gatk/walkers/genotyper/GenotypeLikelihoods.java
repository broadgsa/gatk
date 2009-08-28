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
    //public abstract int add(char ref, char read, byte qual);
}
