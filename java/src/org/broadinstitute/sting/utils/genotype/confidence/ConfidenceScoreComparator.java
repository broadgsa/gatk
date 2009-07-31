package org.broadinstitute.sting.utils.genotype.confidence;

import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.Comparator;


/**
 * 
 * @author aaron 
 * 
 * Class ConfidenceScoreComparator
 *
 * A descriptions should go here. Blame aaron if it's missing.
 */

public class ConfidenceScoreComparator implements Comparator<Genotype> {
    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        return genotype.getConfidenceScore().compareTo(genotype1.getConfidenceScore());
    }
}