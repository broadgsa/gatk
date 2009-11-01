package org.broadinstitute.sting.utils.genotype.vcf;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFGenotypeCall
 *         <p/>
 *         The implementation of the genotype interface, specific to VCF
 */
public class VCFGenotypeCall implements Genotype, ReadBacked, PosteriorsBacked, SampleBacked {
    private final char mRefBase;
    private final GenomeLoc mLocation;

    private List<SAMRecord> mReads;
    private double[] mPosteriors;


    // the reference genotype, the best genotype, and the next best genotype, lazy loaded
    private DiploidGenotype mRefGenotype = null;
    private DiploidGenotype mBestGenotype = null;
    private DiploidGenotype mNextGenotype = null;

    // the sample name, used to propulate the SampleBacked interface
    private String mSampleName;

    /**
     * Generate a single sample genotype object
     *
     */
    public VCFGenotypeCall(char ref, GenomeLoc loc) {
        mRefBase = ref;
        mLocation = loc;

        // fill in empty data
        mPosteriors = new double[10];
        Arrays.fill(mPosteriors, Double.MIN_VALUE);
        mSampleName = "";
        mReads = new ArrayList<SAMRecord>();
    }

    public VCFGenotypeCall(char ref, GenomeLoc loc, String genotype, double negLog10PError, String sample) {
        mRefBase = ref;
        mLocation = loc;
        mBestGenotype = DiploidGenotype.valueOf(genotype);
        mRefGenotype = DiploidGenotype.createHomGenotype(ref);
        mSampleName = sample;

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

        mReads = new ArrayList<SAMRecord>();
    }

    public void setReads(List<SAMRecord> reads) {
        mReads = new ArrayList<SAMRecord>(reads);
    }

    public void setPosteriors(double[] posteriors) {
        mPosteriors = posteriors;
    }

    public void setSampleName(String name) {
        mSampleName = name;
    }


    @Override
    public boolean equals(Object other) {
        lazyEval();

        if (other == null)
            return false;
        if (other instanceof VCFGenotypeCall) {
            VCFGenotypeCall otherCall = (VCFGenotypeCall) other;

            if (!this.mBestGenotype.equals(otherCall.mBestGenotype))
                return false;
            return (this.mRefBase == otherCall.mRefBase);
        }
        return false;
    }

    public String toString() {
        lazyEval();
        return String.format("%s best=%s cmp=%s ref=%s depth=%d negLog10PError=%.2f",
                             getLocation(), mBestGenotype, mRefGenotype, mRefBase, mReads.size(), getNegLog10PError());
    }

    private void lazyEval() {
        if (mBestGenotype == null) {
            Integer sorted[] = Utils.SortPermutation(mPosteriors);
            mBestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
            mNextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
            mRefGenotype = DiploidGenotype.createHomGenotype(this.getReference());
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
    private DiploidGenotype getBestGenotype() {
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
        return !Utils.dupString(this.getReference(), 2).equals(getBestGenotype().toString());
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
    public Variation toVariation() {
        double bestRef = Math.abs(mPosteriors[getBestGenotype().ordinal()] - mPosteriors[getRefGenotype().ordinal()]);
        return new BasicVariation(this.getBases(), String.valueOf(this.getReference()), 0, this.mLocation, bestRef);
    }

    /**
     * get the SAM records that back this genotype call
     *
     * @return a list of SAM records
     */
    public List<SAMRecord> getReads() {
        return mReads;
    }

    /**
     * get the count of reads
     *
     * @return the number of reads we're backed by
     */
    public int getReadCount() {
        return mReads.size();
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

    /**
     * @return returns the sample name for this genotype
     */
    public String getSampleName() {
        return mSampleName;
    }
}