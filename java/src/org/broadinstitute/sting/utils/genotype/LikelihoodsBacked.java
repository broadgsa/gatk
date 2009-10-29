package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 * Interface LikelihoodsBacked
 *
 * this interface indicates that the genotype is
 * backed up by genotype likelihood information.
 */
public interface LikelihoodsBacked {

    /**
     *   
     * @return the likelihood information for this genotype
     */
    public double[] getLikelihoods();

    /**
     *
     * @param   likelihoods    the likelihoods for this genotype
     */
    public void setLikelihoods(double[] likelihoods);

}
