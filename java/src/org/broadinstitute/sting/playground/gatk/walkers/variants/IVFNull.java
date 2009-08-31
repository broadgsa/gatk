package org.broadinstitute.sting.playground.gatk.walkers.variants;

public class IVFNull implements IndependentVariantFeature {
    public void initialize(String arguments) {}

    public void compute(VariantContextWindow contextWindow) {
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

    public boolean useZeroQualityReads() { return false; }
}
