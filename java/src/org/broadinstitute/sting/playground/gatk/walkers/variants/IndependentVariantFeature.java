package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;

/**
 * Interface for conditionally independent variant features.  
 */
public interface IndependentVariantFeature {
    /**
     * A convenient enumeration for each of the ten genotypes.
     */
    public enum Genotype { AA, AC, AG, AT, CC, CG, CT, GG, GT, TT }

    /**
     * Method so that features can initialize themselves based on a short argument string.
     * At the moment, each feature is responsible for interpreting their own argument string.
     *
     * @param arguments
     */
    public void initialize(String arguments);

    /**
     * Method to compute the result of this feature for each of the ten genotypes.  The return value must
     * be a double array of length 10 (one for each genotype) and the value must be in log10-space.
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return a ten-element array of log-likelihood result of the feature applied to each genotype
     */
    public double[] compute(char ref, LocusContext context);
}