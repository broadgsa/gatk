package org.broadinstitute.sting.utils.genotype;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeLikelihood
 *         <p/>
 *         This class encompasses all the information that is associated with a genotype
 *         and it's likelihood, mainly:
 *         <p/>
 *         Likelihood value
 */
public class Genotype {
    private double mLikelihood = 0.0;
    private double mPosteriorProb = 0.0;
    private String mBases = "";
    private int mPloidy = 2; // assume diploid

    /**
     * construct a genotypeLikelihood, given the bases, the posterior, and the likelihood
     *
     * @param bases      the bases that make up this genotype
     * @param posterior  the posterior probability of this genotype
     * @param likelihood the likelihood of this genotype
     * @param ploidy     the ploidy of this genotype
     */
    public Genotype(String bases, double posterior, double likelihood, int ploidy) {
        this.mPloidy = ploidy;
        if (bases.length() != ploidy) {
            throw new IllegalArgumentException("The number of bases should match the ploidy");
        }
        this.mLikelihood = likelihood;
        this.mBases = bases;
        this.mPosteriorProb = posterior;
    }

    /**
     * construct a genotypeLikelihood, given the bases, the posterior, and the likelihood
     *
     * @param bases      the bases that make up this genotype
     * @param posterior  the posterior probability of this genotype
     * @param likelihood the likelihood of this genotype
     */
    public Genotype(String bases, double posterior, double likelihood) {
        if (bases.length() != mPloidy) {
            throw new IllegalArgumentException("The number of bases should match the ploidy");
        }
        this.mLikelihood = likelihood;
        this.mBases = bases;
        this.mPosteriorProb = posterior;
    }

    /**
     * get the likelihood value
     *
     * @return a double, representing the likelihood
     */
    public double getLikelihood() {
        return mLikelihood;
    }

    /**
     * get the posterior value
     *
     * @return a double, representing the posterior
     */
    public double getPosteriorProb() {
        return mPosteriorProb;
    }

    /**
     * get the bases that represent this
     *
     * @return
     */
    public String getBases() {
        return mBases;
    }

    public int getPloidy() {
        return mPloidy;
    }

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return
     */
    public boolean isHom() {
        if (mBases.length() < 1) throw new UnsupportedOperationException("isHom requires at least one stored base");
        char first = mBases.charAt(0);
        for (char c: mBases.toCharArray()) {
            if (c != first) return false;
        }
        return true;
    }

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return
     */
    public boolean isHet() {
         if (mBases.length() < 1) throw new UnsupportedOperationException("isHom requires at least one stored base");
        char first = mBases.charAt(0);
        for (char c: mBases.toCharArray()) {
            if (c != first) return true;
        }
        return false;
    }

}
