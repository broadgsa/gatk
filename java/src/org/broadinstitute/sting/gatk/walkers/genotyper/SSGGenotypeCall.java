package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.confidence.BayesianConfidenceScore;
import org.broadinstitute.sting.utils.genotype.confidence.ConfidenceScore;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class GenotypeCallImpl
 *         <p/>
 *         The single sample genotypers implementation of the genotype call, which contains
 *          extra information for the various genotype outputs
 */
public class SSGGenotypeCall implements GenotypeCall, GenotypeOutput {
    // our stored genotype locus
    private final String mRefBase;
    private final int mPloidy;
    private TreeMap<Double, Genotype> mGenotypes = new TreeMap();
    private double mLikelihoods[];
    private double bestNext = 0;
    private double bestRef = 0;

    private int readDepth = 0;
    private double rmsMapping = 0;

    /**
     * Generate a single sample genotype object, containing everything we need to represent calls out of a genotyper object
     *
     * @param discoveryMode are we in discovery mode or genotyping mode
     * @param mRefBase      the ref base
     * @param mPloidy       the ploidy
     * @param genotypes     the genotype
     * @param likelihoods   the likelihood
     * @param pileup        the pile-up of reads at the specified locus
     */
    public SSGGenotypeCall(boolean discoveryMode, char mRefBase, int mPloidy, List<Genotype> genotypes, double likelihoods[], ReadBackedPileup pileup) {
        this.mRefBase = String.valueOf(mRefBase).toUpperCase();
        this.mPloidy = mPloidy;
        if (genotypes.size() < 1) throw new IllegalArgumentException("Genotypes list size must be greater than 0");
        String refStr = Utils.dupString(mRefBase, mPloidy).toUpperCase();
        calculateBestNext(discoveryMode, genotypes, likelihoods, refStr);

        this.readDepth = pileup.getReads().size();
        rmsMapping = Math.sqrt(calculateRMS(pileup) / readDepth);
    }


    /**
     * calculate the best and next best, and update the stored confidence scores based on these values
     *
     * @param discoveryMode are we in discovery mode?
     * @param genotypes     the genotypes
     * @param likelihoods   the likelihoods corresponding to the genotypes
     * @param refStr        the reference string, i.e. reb base[ploidy]
     */
    private void calculateBestNext(boolean discoveryMode, List<Genotype> genotypes, double[] likelihoods, String refStr) {
        int index = 0;
        double ref = 0.0;
        double best = Double.NEGATIVE_INFINITY;
        double next = Double.NEGATIVE_INFINITY;
        for (Genotype g : genotypes) {
            if (g.getBases().toUpperCase().equals(refStr)) ref = likelihoods[index];
            if (likelihoods[index] > best) {
                next = best;
                best = likelihoods[index];
            } else if (likelihoods[index] > next) next = likelihoods[index];
            index++;
        }
        bestNext = Math.abs(best - next);
        bestRef = Math.abs(best - ref);
        mLikelihoods = likelihoods;
        index = 0;

        // reset the confidence based on either the discovery mode or the genotype mode
        for (Genotype g : genotypes) {
            double val = (discoveryMode) ? Math.abs(likelihoods[index] - ref) : Math.abs(likelihoods[index] - next);
            ((BasicGenotype) g).setConfidenceScore(new BayesianConfidenceScore(val));
            mGenotypes.put(likelihoods[index], g);
            index++;
        }
    }

    @Override
    public boolean equals(Object other) {
        if(other == null)
            return false;
        if(other instanceof SSGGenotypeCall) {
            SSGGenotypeCall otherCall = (SSGGenotypeCall)other;

            boolean likelihoodsMatch = true;
            for ( int i = 0; i < mLikelihoods.length; i++ ) {
                likelihoodsMatch = likelihoodsMatch && MathUtils.compareDoubles(mLikelihoods[i], otherCall.mLikelihoods[i]) == 0;
            }
            //System.out.printf("likelihoodsMatch = %b%n", likelihoodsMatch);

            return this.mRefBase.equals(otherCall.mRefBase) &&
                   this.mPloidy == otherCall.mPloidy &&
                   MathUtils.compareDoubles(this.bestNext, otherCall.bestNext) == 0 &&
                   MathUtils.compareDoubles(this.bestRef, otherCall.bestRef) == 0 &&
                   this.readDepth == otherCall.readDepth &&
                   MathUtils.compareDoubles(this.rmsMapping, otherCall.rmsMapping) == 0 &&
                   likelihoodsMatch;
        }
        return false;
    }

    public String toString() {
        return String.format("%s ref=%s depth=%d rmsMAPQ=%.2f bestVRef=%.2f bestVNext=%.2f likelihoods=%s",
                getLocation(), mRefBase, readDepth, rmsMapping, bestRef, bestNext, Arrays.toString(mLikelihoods));
    }

    /**
     * calculate the rms , given the read pileup
     *
     * @param pileup
     *
     * @return
     */
    private double calculateRMS(ReadBackedPileup pileup) {
        List<SAMRecord> reads = pileup.getReads();
        int[] qualities = new int[reads.size()];
        for (int i=0; i < reads.size(); i++)
            qualities[i] = reads.get(i).getMappingQuality();
        return MathUtils.rms(qualities);
    }

    /**
     * gets the reference base
     *
     * @return the reference base we represent
     */
    @Override
    public char getReferencebase() {
        return mRefBase.charAt(0);
    }

    /**
     * check to see if this called genotype is a variant, i.e. not homozygous reference
     *
     * @return true if it's not hom ref, false otherwise
     */
    @Override
    public boolean isVariation() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).isVariant(getReferencebase());
    }

    /**
     * get the confidence score
     *
     * @return get the confidence score that we're based on
     */
    @Override
    public ConfidenceScore getConfidenceScore() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).getConfidenceScore();
    }

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    @Override
    public String getBases() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).getBases();
    }

    /**
     * get the ploidy
     *
     * @return the ploidy value
     */
    @Override
    public int getPloidy() {
        return this.mPloidy;
    }

    /**
     * Returns true if both observed alleles are the same (regardless of whether they are ref or alt)
     *
     * @return true if we're homozygous, false otherwise
     */
    @Override
    public boolean isHom() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).isHom();
    }

    /**
     * Returns true if observed alleles differ (regardless of whether they are ref or alt)
     *
     * @return true if we're het, false otherwise
     */
    @Override
    public boolean isHet() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).isHet();
    }

    /**
     * Location of this genotype on the reference (on the forward strand). If the allele is insertion/deletion, the first inserted/deleted base
     * is located right <i>after</i> the specified location
     *
     * @return position on the genome wrapped in GenomeLoc object
     */
    @Override
    public GenomeLoc getLocation() {
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).getLocation();
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
        return mGenotypes.get(mGenotypes.descendingKeySet().first()).isVariant(ref);
    }

    /**
     * return this genotype as a variant
     *
     * @return
     */
    @Override
    public Variant toVariant() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * get the genotypes, sorted in asscending order by their ConfidenceScores (the best
     * to the worst ConfidenceScores)
     *
     * @return a list of the genotypes, sorted by ConfidenceScores
     */
    @Override
    public List<Genotype> getGenotypes() {
        List<Genotype> newList = new ArrayList<Genotype>();
        newList.addAll(this.mGenotypes.values());
        return newList;
    }

    /**
     * get the genotypes sorted lexigraphically
     *
     * @return a list of the genotypes sorted lexi
     */
    @Override
    public List<Genotype> getLexigraphicallySortedGenotypes() {
        List<Genotype> newList = new ArrayList<Genotype>();
        newList.addAll(this.mGenotypes.values());
        Collections.sort(newList, new LexigraphicalComparator());
        return newList;
    }

    /**
     * return the likelihoods as a double array, in lexographic order
     *
     * @return the likelihoods
     */
    public double[] getLikelihoods() {
        return this.mLikelihoods;
    }

    public int getReadDepth() {
        return readDepth;
    }

    public double getRmsMapping() {
        return rmsMapping;
    }

    public double getBestNext() {
        return bestNext;
    }

    public double getBestRef() {
        return bestRef;
    }
}


    