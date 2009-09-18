package org.broadinstitute.sting.utils.genotype;

import java.util.List;

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
     * @return a map in lexigraphical order of the likelihoods
     */
    public List<Genotype> getGenotypes();

    /**
     * get the likelihoods
     *
     * @return a map in lexigraphical order of the likelihoods
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
