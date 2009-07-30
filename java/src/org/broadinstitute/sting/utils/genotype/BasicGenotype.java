package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.confidence.ConfidenceScore;


/**
 * @author aaron
 *         <p/>
 *         Class BasicGenotype
 *         <p/>
 *         A basic implementation of the genotype interface
 */
public class BasicGenotype implements Genotype {

    // the bases that represent this genotype
    private String mBases = "";

    // the ploidy, assume 2 unless told otherwise
    private int mPloidy = 2;

    // the confidence score
    protected ConfidenceScore mConfidenceScore;

    // our location
    private GenomeLoc mLoc;

    /**
     * construct a genotypeLikelihood, given the bases, the confidence score, and the ploidy
     *
     * @param loc    the location of the genotype
     * @param bases  the bases that make up this genotype
     * @param ploidy the ploidy of this genotype
     * @param score  the confidence score
     */
    public BasicGenotype(GenomeLoc loc, String bases, int ploidy, ConfidenceScore score) {
        this.mPloidy = ploidy;
        if (bases.length() != ploidy) {
            throw new IllegalArgumentException("The number of bases should match the ploidy");
        }
        this.mBases = bases;
        this.mConfidenceScore = score;
        this.mLoc = loc;
    }

    /**
     * construct a genotypeLikelihood, given the bases and the confidence score, and assume the
     * ploidy is 2.
     *
     * @param loc   the location of the genotype
     * @param bases the bases that make up this genotype
     * @param score the confidence score
     */
    public BasicGenotype(GenomeLoc loc, String bases, ConfidenceScore score) {
        if (bases.length() != mPloidy) {
            throw new IllegalArgumentException("The number of bases should match the ploidy");
        }
        this.mBases = bases;
        this.mConfidenceScore = score;
        this.mLoc = loc;
    }

    /**
     * get the confidence score
     *
     * @return get the confidence score that we're based on
     */
    public ConfidenceScore getConfidenceScore() {
        return this.mConfidenceScore;
    }

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    public String getBases() {
        return mBases;
    }

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    public int getPloidy() {
        return mPloidy;
    }

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    public boolean isHom() {
        if (mBases.length() < 1) throw new UnsupportedOperationException("isHom requires at least one stored base");
        char first = mBases.charAt(0);
        for (char c : mBases.toCharArray()) {
            if (c != first) return false;
        }
        return true;
    }

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    public boolean isHet() {
        if (mBases.length() < 1) throw new UnsupportedOperationException("isHom requires at least one stored base");
        char first = mBases.charAt(0);
        for (char c : mBases.toCharArray()) {
            if (c != first) return true;
        }
        return false;
    }

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * returns true if the genotype is a point genotype, false if it's a indel / deletion
     *
     * @return true is a SNP
     */
    @Override
    public boolean isPointGenotype() {
        return true;
    }

    /**
     * given the reference, are we a variant? (non-ref)
     *
     * @param ref the reference base
     *
     * @return true if we're a variant
     */
    @Override
    public boolean isVariant(char ref) {
        String ret = Utils.dupString(ref, this.getPloidy());
        return !this.getBases().equals(ret);
    }

    /**
     * return this genotype as a variant
     *
     * @return
     */
    @Override
    public Variant toVariant() {
        return null;
    }

    /**
     * set the confidence score
     * @param confidenceScore
     */
    public void setConfidenceScore(ConfidenceScore confidenceScore) {
        this.mConfidenceScore = confidenceScore;
    }

}
