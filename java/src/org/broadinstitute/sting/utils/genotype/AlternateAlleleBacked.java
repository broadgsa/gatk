package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface AlternateAlleleBacked
 *
 * this interface indicates that the genotype is
 * backed up by alternate allele information.
 */
public interface AlternateAlleleBacked {

    /**
     *
     * @return returns the alternate allele for this genotype
     */
    public char getAlternateAllele();

    /**
     *
     * @param   alt    the alternate allele base for this genotype
     */
    public void setAlternateAllele(char alt);

}