package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 * Interface SampleBacked
 *
 * this interface indicates that the genotype is
 * backed up by sample information.
 */
public interface SampleBacked {

    /**
     *
     * @return returns the sample name for this genotype
     */
    public String getSampleName();

    /**
     *
     * @param   name    the sample name for this genotype
     */
    public void setSampleName(String name);
    
}
