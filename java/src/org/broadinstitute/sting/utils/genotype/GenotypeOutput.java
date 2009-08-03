package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 *         <p/>
 *         Interface GenotypeOutputFormat
 *         <p/>
 *         This interface is in adition to the GenotypeCall interface,
 *         but adds the required functionality that output (GenotypeWriter) interfaces
 *         need.  It's fair game to return 0 or any value for these fields, but the methods
 *         are all used by various output writers
 */
public interface GenotypeOutput extends GenotypeCall {

    /**
     * return the likelihoods as a double array, in lexographic order
     *
     * @return the likelihoods
     */
    public double[] getLikelihoods();

    /**
     * get the depth of the reads at this location
     * @return the depth
     */
    public int getReadDepth();

    /**
     * get the rms mapping qualities
     * @return
     */
    public double getRmsMapping();

    /**
     * get the best to the next best genotype
     * @return
     */
    public double getBestNext();

    /**
     * get the best compaired to the reference score
     * @return
     */
    public double getBestRef();
}
