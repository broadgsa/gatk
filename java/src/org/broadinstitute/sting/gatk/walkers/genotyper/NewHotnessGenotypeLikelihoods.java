package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;
import edu.mit.broad.picard.util.MathUtil;

import java.util.Arrays;
import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
public class NewHotnessGenotypeLikelihoods extends GenotypeLikelihoods implements Cloneable {
    private final static int FIXED_PLOIDY = 2;
    private final static int MAX_PLOIDY = FIXED_PLOIDY + 1;

    //
    // The three fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] likelihoods = null;
    protected double[] posteriors = null;

    private DiploidGenotypePriors priors = null;

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     */
    public NewHotnessGenotypeLikelihoods() {
        this.priors = new DiploidGenotypePriors();
        initialize();
    }

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     */
    public NewHotnessGenotypeLikelihoods(DiploidGenotypePriors priors) {
        this.priors = priors;
        initialize();
    }

    /**
     * Cloning of the object
     * @return
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        NewHotnessGenotypeLikelihoods c = (NewHotnessGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.likelihoods = likelihoods.clone();
        c.posteriors = posteriors.clone();
        return c;
    }

    private void initialize() {
        likelihoods = zeros.clone();            // likelihoods are all zeros
        posteriors = priors.getPriors().clone();            // posteriors are all the priors
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return
     */
    public double[] getLikelihoods() {
        return likelihoods;
    }

    /**
     * Returns the likelihood associated with DiploidGenotype g
     * @param g
     * @return log10 likelihood as a double
     */
    public double getLikelihood(DiploidGenotype g) {
        return getLikelihoods()[g.ordinal()];
    }

    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return raw log10 (not-normalized posteriors) as a double array
     */
    public double[] getPosteriors() {
        return posteriors;
    }

    /**
     * Returns the posterior associated with DiploidGenotype g
     * @param g
     * @return raw log10 (not-normalized posteror) as a double
     */
    public double getPosterior(DiploidGenotype g) {
        return getPosteriors()[g.ordinal()];
    }

    public DiploidGenotypePriors getPriorObject() {
        return priors;
    }

    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors.getPriors();
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    /**
     * Are we ignoring Q0 bases during calculations?
     * @return
     */
    public boolean isFilteringQ0Bases() {
        return filterQ0Bases;
    }

    /**
     * Enable / disable filtering of Q0 bases.  Enabled by default
     *
     * @param filterQ0Bases
     */
    public void filterQ0Bases(boolean filterQ0Bases) {
        this.filterQ0Bases = filterQ0Bases;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // add() -- the heart of 
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase
     * @param qualityScore
     * @return 1 if the base was considered good enough to add to the likelihoods (not Q0 or 'N', for example)
     */
    public int add(char observedBase, byte qualityScore, boolean fwdStrand) {
        return add(observedBase, qualityScore, fwdStrand, true, false);
    }

    public int add(char observedBase, byte qualityScore, boolean fwdStrand, boolean enableCache, boolean verbose) {
        if ( badBase(observedBase) ) {
            throw new RuntimeException(String.format("BUG: unexpected base %c with Q%d passed to GenotypeLikelihoods", observedBase, qualityScore));
        }

        // TODO -- we should probably filter Q1 bases too
        if (qualityScore <= 0) {
            if ( isFilteringQ0Bases() ) {
                return 0;
            } else {
                qualityScore = 1;
            }
        }

        // Handle caching if requested.  Just look up the cached result if its available, or compute and store it
        NewHotnessGenotypeLikelihoods cached = null;
        if ( enableCache ) {
            if ( ! inCache(observedBase, qualityScore, FIXED_PLOIDY, fwdStrand) ) {
                cached = calculateCachedGenotypeLikelihoods(observedBase, qualityScore, FIXED_PLOIDY, fwdStrand);
            } else {
                cached = getCachedGenotypeLikelihoods(observedBase, qualityScore, FIXED_PLOIDY, fwdStrand);
            }
        }

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihoodCalc = ! enableCache ? log10PofObservingBaseGivenGenotype(observedBase, g, qualityScore, fwdStrand) : 0.0;
            double likelihoodCached = enableCache ? cached.likelihoods[g.ordinal()] : 0.0;
            double likelihood = enableCache ? likelihoodCached : likelihoodCalc;

            //if ( enableCache && likelihoodCalc != 0.0 && MathUtils.compareDoubles(likelihoodCached, likelihoodCalc) != 0 ) {
            //    System.out.printf("ERROR: Likelihoods not equal is cache=%f != calc=%f for %c %d %s%n",
            //            likelihoodCached, likelihoodCalc, observedBase, qualityScore, g.toString());
            //}
            if ( verbose )
                System.out.printf("  L(%c | G=%s, Q=%d, S=%s) = %f%n", observedBase, g, qualityScore, fwdStrand ? "+" : "-", likelihood);

            likelihoods[g.ordinal()] += likelihood;
            posteriors[g.ordinal()] += likelihood;
        }

        if ( verbose ) {
            for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
            System.out.println();
            for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", likelihoods[g.ordinal()]); }
            System.out.println();
        }

        return 1;
    }

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean verbose) {
        int n = 0;

        for (int i = 0; i < pileup.getReads().size(); i++) {
            SAMRecord read = pileup.getReads().get(i);
            int offset = pileup.getOffsets().get(i);
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];
            if ( ! ignoreBadBases || ! badBase(base) ) {
                n += add(base, qual, ! read.getReadNegativeStrandFlag(), true, verbose);
            }
        }

        return n;
    }

    public int add(ReadBackedPileup pileup, boolean ignoreBadBases) {
        return add(pileup, ignoreBadBases, false);
    }


    public int add(ReadBackedPileup pileup) {
        return add(pileup, false);
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // caching routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------
    final static NewHotnessGenotypeLikelihoods[][][][] cache = new NewHotnessGenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_QUAL_SCORE][MAX_PLOIDY][2];
    static int cacheSize = 0;

    private int strandIndex(boolean fwdStrand) {
        return fwdStrand ? 0 : 1;
    }

    private boolean inCache( char observedBase, byte qualityScore, int ploidy, boolean fwdStrand ) {
        return cache[BaseUtils.simpleBaseToBaseIndex(observedBase)][qualityScore][ploidy][strandIndex(fwdStrand)] != null;
    }

    private NewHotnessGenotypeLikelihoods getCachedGenotypeLikelihoods( char observedBase, byte qualityScore, int ploidy, boolean fwdStrand ) {
        if ( ! inCache(observedBase, qualityScore, ploidy, fwdStrand ) )
            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base=%c, qual=%d, ploidy=%d, fwdstrand=%b",
                    observedBase, qualityScore, ploidy, fwdStrand));

        return cache[BaseUtils.simpleBaseToBaseIndex(observedBase)][qualityScore][ploidy][strandIndex(fwdStrand)];
    }

    private NewHotnessGenotypeLikelihoods calculateCachedGenotypeLikelihoods( char observedBase, byte qualityScore, int ploidy, boolean fwdStrand ) {
        if ( inCache(observedBase, qualityScore, ploidy, fwdStrand ) )
            throw new RuntimeException(String.format("BUG: trying to set an already set cached genotype likelihood at base=%c, qual=%d, ploidy=%d, fwdstrand=%b",
                    observedBase, qualityScore, ploidy, fwdStrand));

        // create a new genotype likelihoods object and add this single base to it -- now we have a cached result
        try {
            NewHotnessGenotypeLikelihoods g = (NewHotnessGenotypeLikelihoods)this.clone();
            g.initialize();
            g.add(observedBase, qualityScore, fwdStrand, false, false);

            cache[BaseUtils.simpleBaseToBaseIndex(observedBase)][qualityScore][ploidy][strandIndex(fwdStrand)] = g;
            cacheSize++;

            //System.out.printf("Caching %c %d %d %d %b (%d total entries)%n", observedBase, BaseUtils.simpleBaseToBaseIndex(observedBase), qualityScore, ploidy, fwdStrand, cacheSize);
            return g;
        } catch ( CloneNotSupportedException e ) {
            throw new RuntimeException(e);
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // helper routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Returns true when the observedBase is considered bad and shouldn't be processed by this object.  A base
     * is considered bad if:
     *
     *   Criterion 1: observed base isn't a A,C,T,G or lower case equivalent
     *
     * @param observedBase
     * @return true if the base is a bad base
     */
    private boolean badBase(char observedBase) {
        return BaseUtils.simpleBaseToBaseIndex(observedBase) == -1;
    }

    /**
     * Return a string representation of this object in a moderately usable form
     *
     * @return
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for (DiploidGenotype g : DiploidGenotype.values()) {
            s.append(String.format("%s %.10f ", g, likelihoods[g.ordinal()]));
			sum += Math.pow(10,likelihoods[g.ordinal()]);
        }
		s.append(String.format(" %f", sum));
        return s.toString();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Validation routines
    //
    //
    // -----------------------------------------------------------------------------------------------------------------

    public boolean validate() {
        return validate(true);
    }

    public boolean validate(boolean throwException) {
        try {
            priors.validate(throwException);

            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                String bad = null;

                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(likelihoods[i]) || ! MathUtils.isNegativeOrZero(likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", likelihoods[i]);
                } else if ( ! MathUtils.wellFormedDouble(posteriors[i]) || ! MathUtils.isNegativeOrZero(posteriors[i]) ) {
                    bad = String.format("Posterior %f is badly formed", posteriors[i]);
                }

                if ( bad != null ) {
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    //
    // Hearty math calculations follow
    //
    //    -- these should not be messed with unless you know what you are doing
    //
    // -----------------------------------------------------------------------------------------------------------------

    /**
     * Works for any ploidy (in organization).
     *
     * @param observedBase
     * @param g
     * @param qual
     * @return
     */
    protected double log10PofObservingBaseGivenGenotype(char observedBase, DiploidGenotype g, byte qual, boolean fwdStrand) {
        if (qual == 0) { // zero quals are wrong
            throw new RuntimeException(String.format("Unexpected Q0 base discovered in calculateAlleleLikelihood: %c %s %d", observedBase, g, qual));
        }

        // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
        double p_base = 0.0;
        p_base += pow(10, log10PofObservingBaseGivenChromosome(observedBase, g.base1, qual, FIXED_PLOIDY, fwdStrand));
        p_base += pow(10, log10PofObservingBaseGivenChromosome(observedBase, g.base2, qual, FIXED_PLOIDY, fwdStrand));
        return log10(p_base);
    }

    protected double log10PofObservingBaseGivenChromosome(char observedBase, char chromBase, byte qual, int ploidy, boolean fwdStrand) {
        double logP = 0.0;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + log10PofTrueBaseGivenMiscall(observedBase, chromBase, fwdStrand);
        }

        // adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
        //logP -= log10N[ploidy];
        logP -= log10(ploidy);

        //System.out.printf("%c %c %d %d => %f%n", observedBase, chromBase, qual, ploidy, logP);
        return logP;
    }

    /**
     * Simple log10(3) cached value
     */
    protected static final double log103 = log10(3.0);

    protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, boolean fwdStrand) {
        return -log103; // currently equivalent to e/3 model
    }

    //
    // Constant static data
    //
    private final static double[] zeros = new double[DiploidGenotype.values().length];

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            zeros[g.ordinal()] = 0.0;
        }
    }
}