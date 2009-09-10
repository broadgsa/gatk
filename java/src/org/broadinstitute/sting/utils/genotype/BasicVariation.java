package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: aaronmckenna
 * Date: Sep 9, 2009
 * Time: 9:32:34 PM
 *
 * a basic implementation of variant
 */
public class BasicVariation implements Variation {

    protected final String mBases;
    protected final char mRef;
    protected final int mLength;
    protected final GenomeLoc mLocation;
    protected final double mConfidence;
    /**
     * the constructor
     * @param bases the bases that this variant represents
     * @param reference the reference base
     * @param length are we a single base variant, or a indel/deletion? length is negitive for an indel,
     * positive for a indel, and 0 for a substitution
     */
    public BasicVariation(String bases, char reference, int length, GenomeLoc location, double confidence) {
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
        if (mLength >0) return VARIANT_TYPE.INDEL;
        if (mLength <0) return VARIANT_TYPE.DELETION;
        return (isSNP()) ? VARIANT_TYPE.SNP : VARIANT_TYPE.REFERENCE;
    }

    @Override
    public boolean isSNP() {
        if (mLength == 0 && Utils.dupString(mRef,2).equals(mBases)) return true;
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
    public String getBases() {
        return mBases;
    }

    @Override
    public GenomeLoc getLocation() {
        return mLocation;
    }

    @Override
    public String getReference() {
        return String.valueOf(mRef);
    }

    @Override
    public boolean isHet() {
        if (mLength == 0 && mBases.charAt(0) != mBases.charAt(1)) return true;
        return false;
    }

    @Override
    public boolean isHom() {
        return !isHet();
    }

    @Override
    public double getNegLog10PError() {
        return mConfidence;
    }

    @Override
    public boolean isReference() {
        if (mLength == 0 && mBases.charAt(0) == mRef && mRef != mBases.charAt(1)) return true;
        return false;
    }

    @Override
    public char getAlternateBase() {
        if (mLength != 0) {
            throw new UnsupportedOperationException("Unable to get alternate base for indel / deletion");
        }
        if (mBases.charAt(0) == mRef) return mBases.charAt(1);
        else return mBases.charAt(0);
    }
}
