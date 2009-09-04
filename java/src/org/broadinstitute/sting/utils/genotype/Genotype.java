package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * @author aaron
 *         <p/>
 *         Class GenotypeLikelihood
 *         <p/>
 *         This class emcompasses all the basic information about a genotype
 */
public interface Genotype {
    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    public double getNegLog10PError();

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    public String getBases();

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    public int getPloidy();

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    public boolean isHom();

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    public boolean isHet();

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation();

    /**
     * returns true if the genotype is a point genotype, false if it's a indel / deletion
     *
     * @return true is a SNP
     */
    public boolean isPointGenotype();

    /**
     * given the reference, are we a variant? (non-ref)
     *
     * @param ref the reference base or bases
     *
     * @return true if we're a variant
     */
    public boolean isVariant(char ref);

    /**
     * get the reference base.
     * @return a character, representing the reference base
     */
    public char getReference();

    /**
     * return this genotype as a variant
     *
     * @return the variant
     */
    public Variant toVariation();

}
