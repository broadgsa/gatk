package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * User: aaron
 * Date: Sep 9, 2009
 * Time: 9:32:34 PM
 * <p/>
 * a basic implementation of the Variation interface.
 */
public class BasicVariation implements Variation {

    // the bases that make up this variant
    protected final String mBases;

    // the reference base
    protected final String mRef;

    // the length of the event, 0 for a SNP, negitive for deletions, positive for insertions
    protected final int mLength;

    // the location on the genome of the event
    protected final GenomeLoc mLocation;

    // our confidence in this event, and a -(log10(Error))
    protected final double mConfidence;

    /**
     * create a basic variation, given the following parameters:
     *
     * @param bases     the bases that this variant represents
     * @param reference the reference bases
     * @param length    are we a single base variant, or a indel/deletion? length is negitive for an indel,
     *                  positive for a indel, and 0 for a substitution
     */
    public BasicVariation(String bases, String reference, int length, GenomeLoc location, double confidence) {
        mBases = bases;
        mRef = reference;
        if (mRef.length() != 1) throw new StingException("The reference must be a single base");
        mLength = length;
        mLocation = location;
        mConfidence = confidence;
    }

    /**
     * we don't know the minor allele freq. is this implementation
     *
     * @return -1.0.  If the freq is less than zero it means we don't know
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return -1.0;
    }

    /**
     * get the type of variation we are
     *
     * @return VARIANT_TYPE
     */
    @Override
    public VARIANT_TYPE getType() {
        if (mLength > 0) return VARIANT_TYPE.INSERTION;
        if (mLength < 0) return VARIANT_TYPE.DELETION;
        return (isSNP()) ? VARIANT_TYPE.SNP : VARIANT_TYPE.REFERENCE;
    }

    @Override
    public boolean isSNP() {
        return ((mLength == 0) && (new HashSet(getAlternateAlleleList()).size() == 1));
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
    public GenomeLoc getLocation() {
        return mLocation;
    }

    @Override
    public String getReference() {
        return (mRef);
    }

    /**
     * are we bi-allelic? In this case we always
     * count the reference as an allele
     */
    @Override
    public boolean isBiallelic() {
        // put the alternate alleles into a set, there may be duplicates (i.e. hom var)
        Set<String> alleles = new HashSet(getAlternateAlleleList());
        return (alleles.size() == 1); // if the alt list contained one unqiue non-ref base, we're biallelic
    }

    @Override
    public double getNegLog10PError() {
        return mConfidence;
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    @Override
    public List<String> getAlternateAlleleList() {
        List<String> list = new ArrayList<String>();
        for (char c : this.mBases.toCharArray())
            if (c != Utils.stringToChar(mRef))
                list.add(String.valueOf(c));
        return list;
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    @Override
    public List<String> getAlleleList() {
        List<String> list = new ArrayList<String>();
        for (char c : this.mBases.toCharArray())
            list.add(String.valueOf(c));
        return list;
    }

    @Override
    public boolean isReference() {
        if (mLength != 0) return false;
        for (String str : getAlleleList())
            if (!str.equals(mRef)) return false;
        return true;
    }

    /**
     * are we an insertion or a deletion? yes, then return true.
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
        if (!this.isBiallelic()) throw new IllegalStateException("we're not biallelic");
        return Utils.stringToChar((new HashSet<String>(getAlternateAlleleList())).iterator().next());
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (!this.isBiallelic()) throw new IllegalStateException("we're not biallelic");
        return Utils.stringToChar(this.mRef);
    }


}
