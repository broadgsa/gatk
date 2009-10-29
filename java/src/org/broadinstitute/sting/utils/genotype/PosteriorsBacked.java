package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 * Interface PosteriorsBacked
 *
 * this interface indicates that the genotype is
 * backed up by genotype posterior information.
 */
public interface PosteriorsBacked {

    /**
     *
     * @return the likelihood information for this genotype
     */
    public double[] getPosteriors();

    /**
     *
     * @param   posteriors    the posteriors for this genotype
     */
    public void setPosteriors(double[] posteriors);

}

