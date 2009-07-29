package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.Comparator;

/**
 * @author aaron
 *         <p/>
 *         Interface Genotype
 *         <p/>
 *         This interface represents the collection of genotypes at a specific locus
 */
public interface GenotypeLocus {

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    public GenomeLoc getLocation();

    /**
     * get the ploidy at this locus
     *
     * @return an integer representing the genotype ploidy at this location
     */
    public int getPloidy();

    /**
     * get the genotypes, sorted in asscending order by their likelihoods (the best
     * to the worst likelihoods)
     *
     * @return a list of the genotypes, sorted by likelihoods
     */
    public List<Genotype> getGenotypes();

    /**
     * get the genotypes and their posteriors
     *
     * @return a list of the genotypes, sorted by poseriors
     */
    public List<Genotype> getPosteriors();

    /**
     * get the genotypes sorted lexigraphically
     *
     * @return a list of the genotypes sorted lexi
     */
    public List<Genotype> getLexigraphicallySortedGenotypes();

    
    /**
     * get the read depth at this position
     *
     * @return the read depth, -1 if it is unknown
     */
    public int getReadDepth();

    /**
     * add a genotype to the collection
     *
     * @param genotype
     *
     * @throws InvalidGenotypeException
     */
    public void addGenotype(Genotype genotype) throws InvalidGenotypeException;

    /**
     * get the root mean square (RMS) of the mapping qualities
     *
     * @return the RMS, or a value < 0 if it's not available
     */
    public double getRMSMappingQuals();

    /**
     * create a variant, given the reference, and a lod score
     *
     * @param refBase the reference base
     * @param score   the threshold to use to determine if it's a variant or not
     *
     * @return a variant object, or null if no genotypes meet the criteria
     */
    public Variant toGenotypeCall(char refBase, ConfidenceScore score);

}


/**
 * the following are helper Comparator classes for the above sort orders, that may be useful
 * for anyone implementing the GenotypeLocus interface
 */
class PosteriorComparator implements Comparator<Genotype> {
    private final Double EPSILON = 1.0e-15;

    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        double diff = genotype.getPosteriorProb() - genotype1.getPosteriorProb();
        if (Math.abs(diff) < (EPSILON * Math.abs(genotype.getPosteriorProb())))
            return 0;
        else if (diff < 0)
            return 1;
        else
            return -1; // TODO: BACKWARD NOW
    }
}

class LexigraphicalComparator implements Comparator<Genotype> {
    private final Double EPSILON = 1.0e-15;

    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        return genotype.getBases().compareTo(genotype1.getBases());
    }
}

class LikelihoodComparator implements Comparator<Genotype> {
    private final Double EPSILON = 1.0e-15;

    @Override
    public int compare(Genotype genotype, Genotype genotype1) {
        double diff = genotype.getLikelihood() - genotype1.getLikelihood();
        if (Math.abs(diff) < (EPSILON * Math.abs(genotype.getLikelihood())))
            return 0;
        else if (diff < 0)
            return 1;  // TODO: BACKWARD NOW
        else
            return -1;
    }
}