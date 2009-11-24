package org.broadinstitute.sting.utils.genotype;


/**
 * @author ebanks
 *         <p/>
 *         Class GenotypeLocusData
 *         <p/>
 *         represents the locus specific data associated with a genotype object.
 */
public interface GenotypeLocusData extends Variation {

    /**
     *
     * @param   alt    the alternate allele base for this genotype
     */
    public void addAlternateAllele(String alt);

    /**
     *
     * @param   frequency    the allele frequency for this genotype
     */
    public void setAlleleFrequency(double frequency);

}