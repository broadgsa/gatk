package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Comparator;
import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Interface GenotypeCall
 *         <p/>
 *         Genotype call interface, for indicating that a genotype is
 *         also a genotype call.
 */
public interface GenotypeCall extends Genotype {
    /**
     * gets the reference base
     *
     * @return the reference base we represent
     */
    public char getReferencebase();

    /**
     * check to see if this called genotype is a variant, i.e. not homozygous reference
     * @return true if it's not hom ref, false otherwise
     */
    public boolean isVariation();

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    public GenomeLoc getLocation();

    /**
     * get the genotypes, sorted in asscending order by their ConfidenceScores (the best
     * to the worst ConfidenceScores)
     *
     * @return a list of the genotypes, sorted by ConfidenceScores
     */
    public List<Genotype> getGenotypes();

    /**
     * get the genotypes sorted lexigraphically
     *
     * @return a list of the genotypes sorted lexi
     */
    public List<Genotype> getLexigraphicallySortedGenotypes();

}

class LexigraphicalComparator implements Comparator<Genotype> {
    private final Double EPSILON = 1.0e-15;

    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        return genotype.getBases().compareTo(genotype1.getBases());
    }
}

class ConfidenceScoreSort implements Comparator<Genotype> {
    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        return genotype.getConfidenceScore().compareTo(genotype1.getConfidenceScore());
    }
}