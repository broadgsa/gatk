package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;

public class IVFNull implements IndependentVariantFeature {
    public String getFeatureName() { return "null"; }

    public double[] compute(char ref, LocusContext context) {
        return new double[0];
    }
}
