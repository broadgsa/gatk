package org.broadinstitute.sting.utils.genotype.geli;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.Arrays;


/**
 * @author ebanks
 *         <p/>
 *         Class GeliGenotypeCall
 *         <p/>
 *         The implementation of the genotype interface, specific to Geli
 */
public class GeliGenotypeCall extends AlleleConstrainedGenotype implements GenotypeCall, ReadBacked, PosteriorsBacked {
    private final char mRefBase;
    private final GenomeLoc mLocation;

    private ReadBackedPileup mPileup = null;
    private double[] mPosteriors;

    private Variation mVariation = null;

    // the reference genotype, the best genotype, and the next best genotype, lazy loaded
    private DiploidGenotype mRefGenotype = null;
    private DiploidGenotype mBestGenotype = null;
    private DiploidGenotype mNextGenotype = null;

    /**
     * Generate a single sample genotype object
     *
     * @param ref  the ref character
     * @param loc  the genome loc
     */
    public GeliGenotypeCall(char ref, GenomeLoc loc) {
        super(ref);
        mRefBase = ref;
        mLocation = loc;

        // fill in empty data
        mPosteriors = new double[10];
        Arrays.fill(mPosteriors, Double.MIN_VALUE);
    }

    public GeliGenotypeCall(char ref, GenomeLoc loc, String genotype, double negLog10PError) {
        super(ref);
        mRefBase = ref;
        mLocation = loc;
        mBestGenotype = DiploidGenotype.valueOf(genotype);
        mRefGenotype = DiploidGenotype.createHomGenotype(ref);
        mNextGenotype = mRefGenotype;

        // set general posteriors to min double value
        mPosteriors = new double[10];
        Arrays.fill(mPosteriors, Double.MIN_VALUE);

        // set the ref to PError
        mPosteriors[mRefGenotype.ordinal()] = -1.0 * negLog10PError;

        // set the best genotype to zero (need to do this after ref in case ref=best)
        mPosteriors[mBestGenotype.ordinal()] = 0.0;

        // choose a smart next best genotype and set it to PError
        if ( mBestGenotype == mRefGenotype )
            mNextGenotype = DiploidGenotype.valueOf(BaseUtils.simpleComplement(genotype));
        else
            mNextGenotype = mRefGenotype;
        mPosteriors[mNextGenotype.ordinal()] = -1.0 * negLog10PError;
    }

    public void setPileup(ReadBackedPileup pileup) {
        mPileup = pileup;
    }

    public void setPosteriors(double[] posteriors) {
        mPosteriors = posteriors;
    }

    public void setVariation(Variation variation) {
        this.mVariation = variation;
    }

    public void setGenotype(DiploidGenotype genotype) {
        ; // do nothing: geli uses diploid posteriors to calculate the genotype
    }

    public void setNegLog10PError(double value) {
        ; // do nothing: geli uses diploid posteriors to calculate the P(error)
    }

    @Override
    public boolean equals(Object other) {
        lazyEval();

        if (other == null)
            return false;
        if (other instanceof GeliGenotypeCall) {
            GeliGenotypeCall otherCall = (GeliGenotypeCall) other;

            if (!this.mBestGenotype.equals(otherCall.mBestGenotype))
                return false;
            return (this.mRefBase == otherCall.mRefBase);
        }
        return false;
    }

    public String toString() {
        lazyEval();
        return String.format("%s best=%s cmp=%s ref=%s depth=%d negLog10PError=%.2f",
                             getLocation(), mBestGenotype, mRefGenotype, mRefBase, getReadCount(), getNegLog10PError());
    }

    private void lazyEval() {
        if (mBestGenotype == null) {
            char ref = this.getReference();
            char alt = this.getAlternateAllele();

            mRefGenotype = DiploidGenotype.createHomGenotype(ref);

            // if we are constrained to a single alternate allele, use only that one
            if ( alt != AlleleConstrainedGenotype.NO_CONSTRAINT ) {
                DiploidGenotype hetGenotype = ref < alt ? DiploidGenotype.valueOf(String.valueOf(ref) + String.valueOf(alt)) : DiploidGenotype.valueOf(String.valueOf(alt) + String.valueOf(ref));
                DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(alt);
                boolean hetOverHom = mPosteriors[hetGenotype.ordinal()] > mPosteriors[homGenotype.ordinal()];
                boolean hetOverRef = mPosteriors[hetGenotype.ordinal()] > mPosteriors[mRefGenotype.ordinal()];
                boolean homOverRef = mPosteriors[homGenotype.ordinal()] > mPosteriors[mRefGenotype.ordinal()];
                if ( hetOverHom ) {
                    mBestGenotype = (hetOverRef ? hetGenotype : mRefGenotype);
                    mNextGenotype = (!hetOverRef ? hetGenotype : (homOverRef ? homGenotype : mRefGenotype));
                } else {
                    mBestGenotype = (homOverRef ? homGenotype : mRefGenotype);
                    mNextGenotype = (!homOverRef ? homGenotype : (hetOverRef ? hetGenotype : mRefGenotype));
                }
            } else {
                Integer sorted[] = Utils.SortPermutation(mPosteriors);
                mBestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
                mNextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
            }
        }
    }


    /**
     * get the confidence we have
     *
     * @return get the one minus error value
     */
    public double getNegLog10PError() {
        return Math.abs(mPosteriors[getBestGenotype().ordinal()] - mPosteriors[getNextBest().ordinal()]);
    }

    // get the best genotype
    protected DiploidGenotype getBestGenotype() {
        lazyEval();
        return mBestGenotype;
    }

    // get the alternate genotype
    private DiploidGenotype getNextBest() {
        lazyEval();
        return mNextGenotype;
    }

    // get the ref genotype
    private DiploidGenotype getRefGenotype() {
        lazyEval();
        return mRefGenotype;
    }

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    public String getBases() {
        return getBestGenotype().toString();
    }

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    public int getPloidy() {
        return 2;
    }

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    public boolean isHom() {
        return getBestGenotype().isHom();
    }

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    public boolean isHet() {
        return !isHom();
    }

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    public GenomeLoc getLocation() {
        return this.mLocation;
    }

    /**
     * returns true if the genotype is a point genotype, false if it's a indel / deletion
     *
     * @return true is a SNP
     */
    public boolean isPointGenotype() {
        return true;
    }

    /**
     * given the reference, are we a variant? (non-ref)
     *
     * @param ref the reference base or bases
     *
     * @return true if we're a variant
     */
    public boolean isVariant(char ref) {
        return !Utils.dupString(ref, 2).equals(getBestGenotype().toString());
    }

    /**
     * are we a variant? (non-ref)
     *
     * @return true if we're a variant
     */
    public boolean isVariant() {
        return isVariant(mRefBase);
    }

    /**
     *
     * @return return this genotype as a variant
     */
    public Variation toVariation(char ref) {
        if ( mVariation == null ) {
            BasicGenotypeBackedVariation var = new BasicGenotypeBackedVariation(ref, mLocation, isVariant() ? Variation.VARIANT_TYPE.SNP : Variation.VARIANT_TYPE.REFERENCE);
            double confidence = Math.abs(mPosteriors[getBestGenotype().ordinal()] - mPosteriors[getRefGenotype().ordinal()]);
            var.setConfidence(confidence);
            if ( isVariant() )
                var.addAlternateAllele(Character.toString(mBestGenotype.base1 != ref ? mBestGenotype.base1 : mBestGenotype.base2));
            mVariation = var;
        }
        return mVariation;
    }

    /**
     * get the pileup that backs this genotype call
     *
     * @return a pileup
     */
    public ReadBackedPileup getPileup() {
        return mPileup;
    }

    /**
     * get the count of reads
     *
     * @return the number of reads we're backed by
     */
    public int getReadCount() {
        return (mPileup != null ? mPileup.getReads().size() : 0);
    }

    /**
     * get the reference character
     *
     * @return the reference character
     */
    public char getReference() {
        return mRefBase;
    }

    /**
     * get the posteriors
     *
     * @return the posteriors
     */
    public double[] getPosteriors() {
        return mPosteriors;
    }
}