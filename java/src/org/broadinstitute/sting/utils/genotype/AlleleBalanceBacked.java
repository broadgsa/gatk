package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface AlleleBalanceBacked
 *
 * this interface indicates that the genotype is
 * backed up by allele balance information.
 */
public interface AlleleBalanceBacked {

    /**
     *
     * @return returns the ratio of ref/(ref+alt) bases for hets
     */
    public double getRefRatio();

    /**
     *
     * @param   ratio    the ref/(ref+alt) ratio
     */
    public void setRefRatio(double ratio);

    /**
     *
     * @return returns the ratio of (ref+alt)/total bases for hets
     */
    public double getOnOffRatio();

    /**
     *
     * @param   ratio    the (ref+alt)/total ratio
     */
    public void setOnOffRatio(double ratio);
}