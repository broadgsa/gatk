package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.Arrays;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class SSGenotypeCall
 *         <p/>
 *         The single sample implementation of the genotype interface, which contains
 *         extra information for the various genotype outputs
 */
public class SSGenotypeCall implements Genotype, ReadBacked, GenotypesBacked, LikelihoodsBacked {
    private final char mRefBase;
    private final GenotypeLikelihoods mGenotypeLikelihoods;

    // the next two values are filled in on-demand.  Default to -1 since they can never be negitive
    private final GenomeLoc mLocation;
    private final ReadBackedPileup mPileup;

    // if this is null, we were constructed with the intention that we'd represent the best genotype
    private DiploidGenotype mGenotype = null;

    // which genotype to compare to; if we're in discovery mode it's the ref allele, otherwise it's the next best
    private DiploidGenotype mCompareTo = null;

    // are we best vrs ref or best vrs next - for internal consumption only
    private final boolean mBestVrsRef;

    /**
     * Generate a single sample genotype object, containing everything we need to represent calls out of a genotyper object
     *
     * @param discovery are we representing the best vrs next or best vrs ref
     * @param location   the location we're working with
     * @param refBase    the ref base
     * @param gtlh       the genotype likelihoods object
     * @param pileup     the pile-up of reads at the specified locus
     */
    public SSGenotypeCall(boolean discovery, GenomeLoc location, char refBase, GenotypeLikelihoods gtlh, ReadBackedPileup pileup) {
        mBestVrsRef = discovery;
        mRefBase = String.valueOf(refBase).toUpperCase().charAt(0); // a round about way to make sure the ref base is up-case
        mGenotypeLikelihoods = gtlh;
        mLocation = location;
        mPileup = pileup;
    }

    /**
     * Generate a single sample genotype object, containing everything we need to represent calls out of a genotyper object
     *
     * @param discovery are we representing the best vrs next or best vrs ref
     * @param location   the location we're working with
     * @param refBase    the ref base
     * @param gtlh       the genotype likelihoods object
     * @param pileup     the pile-up of reads at the specified locus
     */
    SSGenotypeCall(boolean discovery, GenomeLoc location, char refBase, GenotypeLikelihoods gtlh, ReadBackedPileup pileup, DiploidGenotype genotype) {
        mBestVrsRef = discovery;
        mRefBase = String.valueOf(refBase).toUpperCase().charAt(0); // a round about way to make sure the ref base is up-case
        mGenotypeLikelihoods = gtlh;
        mLocation = location;
        mGenotype = genotype;
        mPileup = pileup;
    }

    @Override
    public boolean equals(Object other) {
        if (other == null)
            return false;
        if (other instanceof SSGenotypeCall) {
            SSGenotypeCall otherCall = (SSGenotypeCall) other;

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
        return String.format("%s ref=%s depth=%d rmsMAPQ=%.2f likelihoods=%s",
                             getLocation(), mRefBase, mPileup.getReads().size(), Arrays.toString(mGenotypeLikelihoods.getLikelihoods()));
    }




    /**
     * get the confidence we have
     *
     * @return get the one minus error value
     */
    @Override
    public double getNegLog10PError() {
        getBestGenotype();
        getAltGenotype();        
        return Math.abs(mGenotypeLikelihoods.getPosterior(mGenotype) - mGenotypeLikelihoods.getPosterior(mCompareTo));
    }

    /**
     * get the best genotype
     */
    public DiploidGenotype getBestGenotype() {
        if (mGenotype == null) {
            Integer sorted[] = Utils.SortPermutation(mGenotypeLikelihoods.getPosteriors());
            mGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
        }
        return mGenotype;
    }

    /**
     * get the alternate genotype
     */
    public DiploidGenotype getAltGenotype() {
        if (mCompareTo == null) {
            if (this.mBestVrsRef) {
                mCompareTo = DiploidGenotype.valueOf(Utils.dupString(this.getReference(),2));
            } else {
                Integer sorted[] = Utils.SortPermutation(mGenotypeLikelihoods.getPosteriors());
                mCompareTo = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
            }
        }
        return mCompareTo;
    }

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    @Override
    public String getBases() {
        return getBestGenotype().toString();
    }

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    @Override
    public int getPloidy() {
        return 2;
    }

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    @Override
    public boolean isHom() {
        return getBestGenotype().isHom();
    }

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    @Override
    public boolean isHet() {
        return !isHom();
    }

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    @Override
    public GenomeLoc getLocation() {
        return this.mLocation;
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
     * @param ref the reference base or bases
     *
     * @return true if we're a variant
     */
    @Override
    public boolean isVariant(char ref) {
        return !Utils.dupString(this.getReference(),2).equals(getBestGenotype().toString());
    }

    /**
     * return this genotype as a variant
     *
     * @return
     */
    public Variant toVariation() {
        return null;  // the next step is to implement the variant system
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
        return new SSGenotypeCall(mBestVrsRef,mLocation,mRefBase,mGenotypeLikelihoods,mPileup,x);
    }

    /**
     * get the likelihood information for this
     *
     * @return
     */
    public double[] getLikelihoods() {
        // TODO: this is wrong, obviously, but we've kept it for now to stay backward compatible with previous calls
        return this.mGenotypeLikelihoods.getPosteriors();
    }


}


    