package org.broadinstitute.sting.utils.genotype;

import java.util.Comparator;


/**
 * 
 * @author aaron 
 * 
 * Class LexigraphicalComparator
 *
 * A descriptions should go here. Blame aaron if it's missing.
 */

public class LexigraphicalComparator implements Comparator<Genotype> {
    private final Double EPSILON = 1.0e-15;

    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        return genotype.getBases().compareTo(genotype1.getBases());
    }
}
