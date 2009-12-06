package org.broadinstitute.sting.utils.genotype;


import java.util.List;

/**
 * @author ebanks
 *         <p/>
 *         Class VariationCall
 *         <p/>
 *         represents locus specific (non-genotype) data associated with a call plus the ability to set genotypes.
 */
public interface VariationCall extends Variation {

    /**
     *
     * @param   alt    the alternate allele base for this variation
     */
    public void addAlternateAllele(String alt);

    /**
     *
     * @return true if the allele frequency has been set
     */
    public boolean hasNonRefAlleleFrequency();

    /**
     *
     * @param   frequency    the allele frequency for this variation
     */
    public void setNonRefAlleleFrequency(double frequency);

    /**
     *
     * @param   calls    the Genotypes for this variation
     */
    public void setGenotypeCalls(List<Genotype> calls);

}