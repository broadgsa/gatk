package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.Arrays;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeCall
 *         <p/>
 *         The implementation of the genotype interface, which contains
 *         extra information for the various genotype outputs
 */
public class GenotypeCall implements Genotype, ReadBacked, GenotypesBacked, LikelihoodsBacked, PosteriorsBacked, SampleBacked {
    private final char mRefBase;
    private final GenotypeLikelihoods mGenotypeLikelihoods;

    // the next two values are filled in on-demand.  Default to -1 since they can never be negitive
    private final GenomeLoc mLocation;
    private final ReadBackedPileup mPileup;

    // if this is null, we were constructed with the intention that we'd represent the best genotype
    private DiploidGenotype mGenotype = null;

    // the reference genotype and the next best genotype, lazy loaded
    private DiploidGenotype mRefGenotype = null;
    private DiploidGenotype mNextGenotype = null;

    // the sample name, used to propulate the SampleBacked interface
    private String mSampleName;

    /**
     * Generate a single sample genotype object, containing everything we need to represent calls out of a genotyper object
     *
     * @param sampleName the sample name
     * @param location   the location we're working with
     * @param refBase    the ref base
     * @param gtlh       the genotype likelihoods object
     * @param pileup     the pile-up of reads at the specified locus
     */
    public GenotypeCall(String sampleName, GenomeLoc location, char refBase, GenotypeLikelihoods gtlh, ReadBackedPileup pileup) {
        mSampleName = sampleName;
        mRefBase = Character.toUpperCase(refBase);
        mGenotypeLikelihoods = gtlh;
        mLocation = location;
        mPileup = pileup;
    }

    /**
     * Generate a single sample genotype object, containing everything we need to represent calls out of a genotyper object
     *
     * @param sampleName the sample name
     * @param location   the location we're working with
     * @param refBase    the ref base
     * @param gtlh       the genotype likelihoods object
     * @param pileup     the pile-up of reads at the specified locus
     */
    GenotypeCall(String sampleName, GenomeLoc location, char refBase, GenotypeLikelihoods gtlh, ReadBackedPileup pileup, DiploidGenotype genotype) {
        mSampleName = sampleName;
        mRefBase = Character.toUpperCase(refBase);
        mGenotypeLikelihoods = gtlh;
        mLocation = location;
        mGenotype = genotype;
        mPileup = pileup;
    }

    @Override
    public boolean equals(Object other) {
        lazyEval();

        if (other == null)
            return false;
        if (other instanceof GenotypeCall) {
            GenotypeCall otherCall = (GenotypeCall) other;

            if (!this.mGenotypeLikelihoods.equals(otherCall.mGenotypeLikelihoods)) {
                return false;
            }
            if (!this.mGenotype.equals(otherCall.mGenotype))
                return false;
            return (this.mRefBase == otherCall.mRefBase) &&
                    this.mPileup.equals(mPileup);
        }
        return false;
    }

    public String toString() {
        lazyEval();
        return String.format("%s best=%s cmp=%s ref=%s depth=%d negLog10PError = %.2f, likelihoods=%s",
                             getLocation(), mGenotype, mRefGenotype, mRefBase, mPileup.getReads().size(),
                             getNegLog10PError(), Arrays.toString(mGenotypeLikelihoods.getLikelihoods()));
    }

    private void lazyEval() {
        // us
        if (mGenotype == null) {
            Integer sorted[] = Utils.SortPermutation(mGenotypeLikelihoods.getPosteriors());
            mGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
        }

        // our comparison
        if (mRefGenotype == null) {
            mRefGenotype = DiploidGenotype.valueOf(Utils.dupString(this.getReference(), 2));
        }
        if (mNextGenotype == null) {
            Integer sorted[] = Utils.SortPermutation(mGenotypeLikelihoods.getPosteriors());
            mNextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
        }
    }


    /**
     * get the confidence we have
     *
     * @return get the one minus error value
     */
    public double getNegLog10PError() {
        return Math.abs(mGenotypeLikelihoods.getPosterior(getBestGenotype()) - mGenotypeLikelihoods.getPosterior(getNextBest()));
    }

    /** get the best genotype */
    private DiploidGenotype getBestGenotype() {
        lazyEval();
        return mGenotype;
    }

    /** get the alternate genotype */
    private DiploidGenotype getNextBest() {
        lazyEval();
        return mNextGenotype;
    }

    /** get the alternate genotype */
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
     * return this genotype as a variant
     *
     * @return
     */
    public Variation toVariation() {
        double bestRef = Math.abs(mGenotypeLikelihoods.getPosterior(getBestGenotype()) - mGenotypeLikelihoods.getPosterior(getRefGenotype()));
        return new BasicVariation(this.getBases(), String.valueOf(this.getReference()), 0, this.mLocation, bestRef);
    }

    /**
     * return the likelihoods as a double array, in lexographic order
     *
     * @return the likelihoods
     */
    public double[] getProbabilities() {
        return this.mGenotypeLikelihoods.getPosteriors();
    }

    /**
     * get the SAM records that back this genotype call
     *
     * @return a list of SAM records
     */
    public List<SAMRecord> getReads() {
        return this.mPileup.getReads();
    }

    /**
     * get the count of reads
     *
     * @return the number of reads we're backed by
     */
    public int getReadCount() {
        return this.mPileup.getReads().size();
    }

    /**
     * get the reference character
     *
     * @return
     */
    public char getReference() {
        return this.mRefBase;
    }

    /**
     * get the likelihoods
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    public Genotype getGenotype(DiploidGenotype x) {
        return new GenotypeCall(mSampleName, mLocation, mRefBase, mGenotypeLikelihoods, mPileup, x);
    }

    /**
     * get the likelihood information for this
     *
     * @return
     */
    public double[] getLikelihoods() {
        return this.mGenotypeLikelihoods.getLikelihoods();
    }


    /**
     * get the likelihood information for this
     *
     * @return
     */
    public double[] getPosteriors() {
        return this.mGenotypeLikelihoods.getPosteriors();
    }

    /** @return returns the sample name for this genotype */
    public String getSampleName() {
        return this.mSampleName;
    }

    /**
     * get the filtering string for this genotype
     *
     * @return a string, representing the genotyping value
     */
    public String getFilteringValue() {
        return "0";
    }
}


    