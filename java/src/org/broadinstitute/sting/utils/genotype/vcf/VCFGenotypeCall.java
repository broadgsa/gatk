package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;


/**
 * @author ebanks
 *         <p/>
 *         Class VCFGenotypeCall
 *         <p/>
 *         The implementation of the genotype interface, specific to VCF
 */
public class VCFGenotypeCall extends AlleleConstrainedGenotype implements GenotypeCall, ReadBacked, SampleBacked, Cloneable {
    private final String mRef;
    private final GenomeLoc mLocation;

    private ReadBackedPileup mPileup = null;
    private int mCoverage = 0;
    private double mNegLog10PError = -1;

    private Variation mVariation = null;

    // the best genotype
    private DiploidGenotype mGenotype = null;

    // the sample name, used to propulate the SampleBacked interface
    private String mSampleName;

    public VCFGenotypeCall(String ref, GenomeLoc loc) {
        super(ref);
        mRef = ref;
        mLocation = loc;

        // fill in empty data
        mGenotype = DiploidGenotype.createHomGenotype(ref.charAt(0));
        mSampleName = "";
    }

    public VCFGenotypeCall(String ref, GenomeLoc loc, DiploidGenotype genotype, double negLog10PError, int coverage, String sample) {
        super(ref);
        mRef = ref;
        mLocation = loc;
        mGenotype = genotype;
        mNegLog10PError = negLog10PError;
        mCoverage = coverage;
        mSampleName = sample;
    }

    public void setPileup(ReadBackedPileup pileup) {
        mPileup = pileup;
    }

    public void setGenotype(DiploidGenotype genotype) {
        mGenotype = genotype;
    }

    public void setNegLog10PError(double negLog10PError) {
        mNegLog10PError = negLog10PError;
    }

    public void setVariation(Variation variation) {
        this.mVariation = variation;
    }

    public void setSampleName(String name) {
        mSampleName = name;
    }


    @Override
    public boolean equals(Object other) {

        if ( other == null || !(other instanceof VCFGenotypeCall) )
            return false;

        VCFGenotypeCall otherCall = (VCFGenotypeCall) other;

        return mGenotype.equals(otherCall.mGenotype) &&
               mNegLog10PError == otherCall.mNegLog10PError &&
               mLocation.equals(otherCall.mLocation) &&
               mRef.equals(otherCall.mRef);
    }

    public String toString() {
        return String.format("%s best=%s ref=%s depth=%d negLog10PError=%.2f",
                             getLocation(), mGenotype, mRef, getReadCount(), getNegLog10PError());
    }

    /**
     * get the confidence we have
     *
     * @return get the one minus error value
     */
    public double getNegLog10PError() {
        return mNegLog10PError;
    }

    // get the best genotype
    protected DiploidGenotype getBestGenotype() {
        return mGenotype;
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

    // You can't make a 'no call' genotype call
    public boolean isNoCall() { return false; }

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
        return !getBestGenotype().isHomRef(ref);
    }

    /**
     * are we a variant? (non-ref)
     *
     * @return true if we're a variant
     */
    public boolean isVariant() {
        return isVariant(mRef.charAt(0));
    }

    /**
     *
     * @return return this genotype as a variant
     */
    public Variation toVariation(char ref) {
        if ( mVariation == null ) {
            VCFVariationCall var = new VCFVariationCall(ref, mLocation, isVariant() ? Variation.VARIANT_TYPE.SNP : Variation.VARIANT_TYPE.REFERENCE);
            var.setConfidence(10 * mNegLog10PError);
            if ( !mGenotype.isHomRef(ref) ) {
                if ( mGenotype.base1 != ref )
                    var.addAlternateAllele(Character.toString(mGenotype.base1));
                if ( mGenotype.isHet() && mGenotype.base2 != ref )
                    var.addAlternateAllele(Character.toString(mGenotype.base2));
            }
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
       return (mCoverage > 0 ? mCoverage : (mPileup != null ? mPileup.size() : 0));
    }

    /**
     * get the reference string
     *
     * @return the reference string
     */
    public String getReference() {
        return mRef;
    }

    /**
     * @return returns the sample name for this genotype
     */
    public String getSampleName() {
        return mSampleName;
    }
}