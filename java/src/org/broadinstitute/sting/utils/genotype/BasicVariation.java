package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Sep 9, 2009
 * Time: 9:32:34 PM
 * <p/>
 * a basic implementation of variant
 */
public class BasicVariation implements Variation {

    protected final String mBases;
    protected final String mRef;
    protected final int mLength;
    protected final GenomeLoc mLocation;
    protected final double mConfidence;

    /**
     * the constructor
     *
     * @param bases     the bases that this variant represents
     * @param reference the reference bases
     * @param length    are we a single base variant, or a indel/deletion? length is negitive for an indel,
     *                  positive for a indel, and 0 for a substitution
     */
    public BasicVariation(String bases, String reference, int length, GenomeLoc location, double confidence) {
        mBases = bases;
        mRef = reference;
        mLength = length;
        mLocation = location;
        mConfidence = confidence;
    }

    @Override
    public double getNonRefAlleleFrequency() {
        return -1.0;
    }

    @Override
    public VARIANT_TYPE getType() {
        if (mLength > 0) return VARIANT_TYPE.INDEL;
        if (mLength < 0) return VARIANT_TYPE.DELETION;
        return (isSNP()) ? VARIANT_TYPE.SNP : VARIANT_TYPE.REFERENCE;
    }

    @Override
    public boolean isSNP() {
        if (mLength == 0) return true;
        return false;
    }

    @Override
    public boolean isInsertion() {
        return (mLength > 0);
    }

    @Override
    public boolean isDeletion() {
        return (mLength < 0);
    }

    @Override
    public String getAlternateBases() {
        return mBases;
    }

    @Override
    public GenomeLoc getLocation() {
        return mLocation;
    }

    @Override
    public String getReference() {
        return (mRef);
    }

    @Override
    public double getNegLog10PError() {
        return mConfidence;
    }

    @Override
    public boolean isReference() {
        if (mLength != 0) return true;
        int refIndex = 0;
        for (char c : mBases.toCharArray()) {
            if (mRef.charAt(refIndex) != c) return false;
        }
        return true;
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false it is.
     *
     * @return true if we're an insertion or deletion
     */
    @Override
    public boolean isIndel() {
        return (isDeletion() || isInsertion());
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (getAlternateBases().charAt(0) == this.getReference().charAt(0))
            return getAlternateBases().charAt(1);
        return getAlternateBases().charAt(0);
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (getAlternateBases().charAt(0) == this.getReference().charAt(0))
            return getAlternateBases().charAt(0);
        return getAlternateBases().charAt(1);
    }


}
