package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;

public interface IndependentVariantFeature {
    public enum Genotype { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    public void initialize(String arguments);
    
    public double[] compute(char ref, LocusContext context);
}
