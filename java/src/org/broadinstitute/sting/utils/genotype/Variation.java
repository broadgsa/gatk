package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

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
        SNP, INDEL, DELETION, REFERENCE // though reference is not really a variant
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    public double getNonRefAlleleFrequency();

    /** @return the VARIANT_TYPE of the current variant */
    public VARIANT_TYPE getType();

    /**
     * are we a SNP? If not we're a Indel/deletion
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
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isDeletion();

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
    public char getReference();

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    public double getNegLog10PError();

    /**
     * are we truely a variant, given a reference
     * @return false if we're a variant(indel, delete, SNP, etc), true if we're not
     */
    public boolean isReference();

    /**
     * gets the alternate bases.  If this is homref, throws an UnsupportedOperationException
     * @return
     */
    public String getAlternateBases();

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false it is.
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel();

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of 
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP();

}
