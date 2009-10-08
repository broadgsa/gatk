package org.broadinstitute.sting.utils.genotype;

/**
 * @author aaron
 *         <p/>
 *         Interface SampleBacked
 *         <p/>
 *         A descriptions should go here. Blame aaron if it's missing.
 */
public interface SampleBacked {

    /**
     *
     * @return returns the sample name for this genotype
     */
    public String getSampleName();

    /**
     * get the filtering string for this genotype
     * @return a string, representing the genotyping value
     */
    public String getFilteringValue();
}
