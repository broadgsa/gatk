package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;

/**
 * @author aaron
 *
 * Interface LikelihoodsBacked
 *
 * indicate that this genotype was called, and contains the other genotypes with their
 * associtated information.
 */
public interface GenotypesBacked {
    /**
     * get the likelihoods
     * @return an array in lexigraphical order of the likelihoods
     */
    public Genotype getGenotype(DiploidGenotype x);
}
