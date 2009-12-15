package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 * Interface PosteriorsBacked
 *
 * This interface indicates that the genotype is backed up by *diploid* genotype posterior information.
 * The posteriors array is expected to have 10 entries (one for each of the possible diploid genotypes).
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

