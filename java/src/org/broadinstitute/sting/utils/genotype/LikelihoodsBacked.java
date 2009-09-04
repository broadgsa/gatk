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
     * get the likelihood information for this 
     * @return
     */
    public double[] getLikelihoods();
}
