package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

/**
 * @author aaron
 *         <p/>
 *         Interface VariantBackedByGenotype
 *         <p/>
 *         this variant is backed by genotypic information
 */
public interface VariantBackedByGenotype {
    /**
     * get the likelihoods
     *
     * @param x the genotype
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    public Genotype getGenotype(DiploidGenotype x);


    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    public boolean hasGenotype(DiploidGenotype x);
}
