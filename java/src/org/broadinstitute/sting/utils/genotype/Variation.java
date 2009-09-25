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
        SNP, INDEL, REFERENCE // though reference is not really a variant
    }

    /**
     * get the frequency of this variant, if we're a variant.  If we're reference this method
     * should return 0.
     *
     * @return double with the stored frequency
     */
    public double getNonRefAlleleFrequency();

    /**
     * A convenience method, for switching over the variation type
     * @return the VARIANT_TYPE of the current variant
     **/
    public VARIANT_TYPE getType();

    /**
     * are we a SNP? If not we're a Indel/deletion or the reference
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
     * get the location that this Variant represents
     *
     * @return a GenomeLoc
     */
    public GenomeLoc getLocation();

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    public String getReference();

    /** are we bi-allelic? */
    public boolean isBiallelic();

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    public double getNegLog10PError();

    /**
     * gets the alternate base.  Use this method if we're biallelic
     *
     * @return
     */
    public String getAlternateBases();

    /**
     * gets the alternate bases.  Use this method if the allele count is greater then 2 (not biallelic)
     *
     * @return
     */
    public List<String> getAlternateBaseList();

    /**
     * are we an insertion or a deletion? yes, then return true.  No? false.
     *
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel();

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
     * @return a char, representing the alternate base
     */
    public char getReferenceForSNP();

}
