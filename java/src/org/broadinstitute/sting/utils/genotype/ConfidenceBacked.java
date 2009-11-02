package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface ConfidenceBacked
 *
 * this interface indicates that the genotype is
 * backed up by sample information.
 */
public interface ConfidenceBacked {

    /**
     *
     * @return returns the confidence for this genotype
     */
    public double getConfidence();

    /**
     *
     * @param   confidence    the confidence for this genotype
     */
    public void setConfidence(double confidence);

}