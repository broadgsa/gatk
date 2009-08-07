package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import edu.mit.broad.picard.util.MathUtil;

import java.util.Arrays;

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
public class NewHotnessGenotypeLikelihoods extends GenotypeLikelihoods {
    private int coverage = 0;

    //
    // The three fundamental data arrays associated with a Genotype Likelhoods object
    //
    private double[] likelihoods = null;
    private double[] posteriors = null;

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

    public int getCoverage() {
        return coverage;
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

    /**
     * Remove call -- for historical purposes only
     */
    @Deprecated
    public int add(char ref, char read, byte qual) {
        return add(read, qual);
    }

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase
     * @param qualityScore
     * @return 1 if the base was considered good enough to add to the likelihoods (not Q0 or 'N', for example)
     */
    public int add(char observedBase, byte qualityScore)
	{
        if ( badBase(observedBase) ) {
            throw new RuntimeException(String.format("BUG: unexpected base %c with Q%d passed to GenotypeLikelihoods", observedBase, qualityScore));
        }

        if (qualityScore <= 0) {
            if ( isFilteringQ0Bases() ) {
                return 0;
            } else {
                qualityScore = 1;
            }
        }

        coverage++;
        for (DiploidGenotype g : DiploidGenotype.values() ) {
			double likelihood = calculateBaseLikelihood(observedBase, g, qualityScore);
            //System.out.printf("Likelihood is %f for %c %d %s%n", likelihood, read, qual, g.toString());
            likelihoods[g.ordinal()] += likelihood;
            posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }

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
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases) {
        int n = 0;

        for (int i = 0; i < pileup.getReads().size(); i++) {
            SAMRecord read = pileup.getReads().get(i);
            int offset = pileup.getOffsets().get(i);
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];
            if ( ! ignoreBadBases || ! badBase(base) ) {
                n += add(base, qual);
            }
        }

        return n;
    }

    public int add(ReadBackedPileup pileup) {
        return add(pileup, false);
    }

    private double calculateBaseLikelihood(char observedBase, DiploidGenotype g, byte qual) {
        if (qual == 0) { // zero quals are wrong
            throw new RuntimeException(String.format("Unexpected Q0 base discovered in calculateAlleleLikelihood: %c %s %d", observedBase, g, qual));
        }

        double p_base = 0.0;

        if ( g.isHomRef(observedBase) ) {
            // hom
            p_base = getOneMinusQual(qual);
        } else if ( g.isHetRef(observedBase) ) {
            // het
            p_base = getOneHalfMinusQual(qual);
        } else {
            // error
            //System.out.printf("%s %b %f %f%n", genotype, h1 != h2, log10Of2_3, log10Of1_3 );
            p_base = qual / -10.0 + log10Of1_3;
        }

        return p_base;
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

    private static double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData3Base[qual];
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