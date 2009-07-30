package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.ReadBackedPileup;


/**
 * 
 * @author aaron 
 * 
 * Class GenotypeFactory
 *
 * Create genotypes, given certain pileup information
 */
public interface GenotypeGenerator {

    /**
     * given all the data associated with a locus, make a genotypeLocus object containing the likelihoods and posterior probs
     * 
     * @param tracker contains the reference meta data for this location, which may contain relevent information like dpSNP or hapmap information
     * @param ref the reference base
     * @param pileup a pileup of the reads, containing the reads and their offsets
     * @return a GenotypeLocus, containing each of the genotypes and their associated likelihood and posterior prob values
     */
    public GenotypeCall callGenotypes(RefMetaDataTracker tracker, char ref, ReadBackedPileup pileup);
}
