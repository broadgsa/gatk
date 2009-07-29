package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;

public class IVFNull implements IndependentVariantFeature {
    public void initialize(String arguments) {}

    public void compute(char ref, LocusContext context) {
    }

    public double[] getLikelihoods() {
        return new double[10];
    }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }
}
