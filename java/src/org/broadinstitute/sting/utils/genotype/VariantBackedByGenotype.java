package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;

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
     * @return an array in lexigraphical order of the likelihoods
     */
    public Genotype getGenotype(DiploidGenotype x);
}
