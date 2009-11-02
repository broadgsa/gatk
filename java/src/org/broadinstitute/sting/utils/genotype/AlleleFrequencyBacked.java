package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface AlleleFrequencyBacked
 *
 * this interface indicates that the genotype is
 * backed up by allele frequency information.
 */
public interface AlleleFrequencyBacked {

    /**
     *
     * @return returns the allele frequency for this genotype
     */
    public double getAlleleFrequency();

    /**
     *
     * @param   frequency    the allele frequency for this genotype
     */
    public void setAlleleFrequency(double frequency);

}