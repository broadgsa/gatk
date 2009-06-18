package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;

public class IVFNull implements IndependentVariantFeature {
    public void initialize(String arguments) {}

    public double[] compute(char ref, LocusContext context) {
        return new double[10];
    }
}
