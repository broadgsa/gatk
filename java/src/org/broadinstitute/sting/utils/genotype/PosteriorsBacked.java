package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 *         <p/>
 *         Interface PosteriorsBacked
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public interface PosteriorsBacked {

    /**
     * get the likelihood information for this
     * @return
     */
    public double[] getPosteriors();
}

