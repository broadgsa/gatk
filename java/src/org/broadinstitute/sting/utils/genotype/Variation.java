package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * @author aaron
 *         <p/>
 *         Interface Variant
 *         <p/>
 *         This class represents a variant
 */
public interface Variation {
    // the types of variants we currently allow
    public enum VARIANT_TYPE {
        SNP, INSERTION, DELETION, REFERENCE // though reference is not really a variant, we need to represent it
    }

    /**
     * @return true if we are bi-allelic?
     */
    public boolean isBiallelic();

    /**
     * get the frequency of this variant, if we're a variant.  If we're reference this method
     * should return 0.  If we can't provide an alternate allele frequency, this should also
     * return 0.
     *
     * WARNING: This method is only valid for biAllelic data, the contract is to check isBiallelic()
     * before calling this method
     *
     * @return double the minor allele frequency
     */
    public double getNonRefAlleleFrequency();

    /**
     * A convenience method, for switching over the variation type
     * @return the VARIANT_TYPE of the current variant
     **/
    public VARIANT_TYPE getType();

    /**
     * are we a SNP? If not we're a Indel/deletion or the reference.  This method must be called  before you use
     * the convenience methods getAlternativeBaseForSNP or getReferenceForSNP, to ensure that you're working with a SNP
     *
     * @return true if we're a SNP
     */
    public boolean isSNP();

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isInsertion();

    /**
     * are we an deletion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isDeletion();

    /**
     * are we a variant that represents the reference allele?
     *
     * @return false if we're a variant(indel, delete, SNP, etc), true if we're hom ref
     */
    public boolean isReference();

    /**
     * are we an insertion or a deletion? yes, then return true.  No? false.
     *
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel();

    /**
     * get the location of this Variant
     *
     * @return a GenomeLoc
     */
    public GenomeLoc getLocation();

    /**
     * get the reference base(s) for this Variant
     *
     * @return the reference base or bases, as a string
     */
    public String getReference();

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    public double getNegLog10PError();


    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    public List<String> getAlternateAlleleList();

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles. If the reference base is not an allele in this varation
     * it will not be in the list (i.e. there is no guarantee that the reference base is in the list).
     *
     * @return an alternate allele list
     */
    public List<String> getAlleleList();

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP();

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the reference base
     */
    public char getReferenceForSNP();

}
