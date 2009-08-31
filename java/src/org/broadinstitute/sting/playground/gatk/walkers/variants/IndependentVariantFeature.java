package org.broadinstitute.sting.playground.gatk.walkers.variants;

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
     * @param arguments the arguments!
     */
    public void initialize(String arguments);

    /**
     * Method to compute the result of this feature for each of the ten genotypes.  The return value must
     * be a double array of length 10 (one for each genotype) and the value must be in log10-space.
     * 
     * @param contextWindow  the context for the given locus
     */
    public void compute(VariantContextWindow contextWindow);

    public double[] getLikelihoods();

    public String getStudyHeader();

    public String getStudyInfo();
}