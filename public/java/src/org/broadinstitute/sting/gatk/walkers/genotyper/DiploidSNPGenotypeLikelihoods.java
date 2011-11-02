/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.fragments.FragmentCollection;
import org.broadinstitute.sting.utils.fragments.FragmentUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.List;

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
public class DiploidSNPGenotypeLikelihoods implements Cloneable {
    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;

    protected final static int FIXED_PLOIDY = 2;
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;
    protected final static double ploidyAdjustment = log10(FIXED_PLOIDY);
    protected final static double log10_3 = log10(3.0);

    protected boolean VERBOSE = false;

    //
    // The fundamental data arrays associated with a Genotype Likelhoods object
    //
    protected double[] log10Likelihoods = null;
    protected double[] log10Posteriors = null;

    protected DiploidSNPGenotypePriors priors = null;

    // TODO: don't calculate this each time through
    protected double log10_PCR_error_3;
    protected double log10_1_minus_PCR_error;

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     */
    public DiploidSNPGenotypeLikelihoods() {
        this.priors = new DiploidSNPGenotypePriors();
        log10_PCR_error_3 = log10(DEFAULT_PCR_ERROR_RATE) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - DEFAULT_PCR_ERROR_RATE);
        setToZero();
    }

    /**
     * Create a new GenotypeLikelhoods object with given priors and PCR error rate for each diploid genotype
     *
     * @param priors          priors
     * @param PCR_error_rate  the PCR error rate
     */
    public DiploidSNPGenotypeLikelihoods(DiploidSNPGenotypePriors priors, double PCR_error_rate) {
        this.priors = priors;
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
        setToZero();
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        DiploidSNPGenotypeLikelihoods c = (DiploidSNPGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        return c;
    }

    protected void setToZero() {
        log10Likelihoods = genotypeZeros.clone();                 // likelihoods are all zeros
        log10Posteriors = priors.getPriors().clone();     // posteriors are all the priors
    }

    /**
     * Returns an array of log10 likelihoods for each genotype, indexed by DiploidGenotype.ordinal values()
     * @return likelihoods array
     */
    public double[] getLikelihoods() {
        return log10Likelihoods;
    }

    /**
     * Returns the likelihood associated with DiploidGenotype g
     * @param g genotype
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
        return log10Posteriors;
    }

    /**
     * Returns the posterior associated with DiploidGenotype g
     * @param g genotpe
     * @return raw log10 (not-normalized posteror) as a double
     */
    public double getPosterior(DiploidGenotype g) {
        return getPosteriors()[g.ordinal()];
    }


    /**
     * Returns an array of posteriors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return normalized posterors as a double array
     */
    public double[] getNormalizedPosteriors() {
        double[] normalized = new double[log10Posteriors.length];
        double sum = 0.0;

        // for precision purposes, we need to add (or really subtract, since everything is negative)
        // the largest posterior value from all entries so that numbers don't get too small
        double maxValue = log10Posteriors[0];
        for (int i = 1; i < log10Posteriors.length; i++) {
            if ( maxValue < log10Posteriors[i] )
                maxValue = log10Posteriors[i];
        }

        // collect the posteriors
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double posterior = Math.pow(10, getPosterior(g) - maxValue);
            normalized[g.ordinal()] = posterior;
            sum += posterior;
        }

        // normalize
        for (int i = 0; i < normalized.length; i++)
            normalized[i] /= sum;

        return normalized;
    }



    public DiploidSNPGenotypePriors getPriorObject() {
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
     * Sets the priors
     * @param priors priors
     */
    public void setPriors(DiploidSNPGenotypePriors priors) {
        this.priors = priors;
        log10Posteriors = genotypeZeros.clone();
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Posteriors[i] = priors.getPriors()[i] + log10Likelihoods[i];
        }
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g genotype
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    // -------------------------------------------------------------------------------------
    //
    // add() routines.  These are the workhorse routines for calculating the overall genotype
    // likelihoods given observed bases and reads.  Includes high-level operators all the
    // way down to single base and qual functions.
    //
    // -------------------------------------------------------------------------------------

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @param minBaseQual               the minimum base quality at which to consider a base valid
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        int n = 0;

        // for each fragment, add to the likelihoods
        FragmentCollection<PileupElement> fpile = pileup.toFragments();

        for ( PileupElement p : fpile.getSingletonReads() )
            n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        for ( List<PileupElement> overlappingPair : fpile.getOverlappingPairs() )
            n += add(overlappingPair, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        return n;
    }

    public int add(PileupElement elt, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        byte obsBase = elt.getBase();

        if ( elt.isReducedRead() ) {
            // reduced read representation
            byte qual = elt.getQual();
            if ( BaseUtils.isRegularBase( elt.getBase() )) {
                add(obsBase, qual, (byte)0, (byte)0, elt.getRepresentativeCount()); // fast calculation of n identical likelihoods
                return elt.getRepresentativeCount(); // we added nObs bases here
            } else // odd bases or deletions => don't use them
                return 0;
        } else {
            byte qual = qualToUse(elt, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
            return qual > 0 ? add(obsBase, qual, (byte)0, (byte)0, 1) : 0;
        }
    }

    public int add(List<PileupElement> overlappingPair, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        final PileupElement p1 = overlappingPair.get(0);
        final PileupElement p2 = overlappingPair.get(1);

        final byte observedBase1 = p1.getBase();
        final byte qualityScore1 = qualToUse(p1, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);
        final byte observedBase2 = p2.getBase();
        final byte qualityScore2 = qualToUse(p2, ignoreBadBases, capBaseQualsAtMappingQual, minBaseQual);

        if ( qualityScore1 == 0 ) {
            if ( qualityScore2 == 0 ) // abort early if we didn't see any good bases
                return 0;
            else {
                return add(observedBase2, qualityScore2, (byte)0, (byte)0);
            }
        } else {
            return add(observedBase1, qualityScore1, observedBase2, qualityScore2);
        }
    }

    /**
     *
     * @param obsBase1
     * @param qual1
     * @param obsBase2
     * @param qual2 can be 0, indicating no second base was observed for this fragment
     * @param nObs The number of times this quad of values was seen.  Generally 1, but reduced reads
     *  can have nObs > 1 for synthetic reads
     * @return
     */
    private int add(byte obsBase1, byte qual1, byte obsBase2, byte qual2, int nObs) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.

        // Just look up the cached result if it's available, or compute and store it
        DiploidSNPGenotypeLikelihoods gl;
        if ( ! inCache(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY) ) {
            gl = calculateCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        } else {
            gl = getCachedGenotypeLikelihoods(obsBase1, qual1, obsBase2, qual2, FIXED_PLOIDY);
        }

        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            log10Likelihoods[g.ordinal()] += likelihood * nObs;
            log10Posteriors[g.ordinal()] += likelihood * nObs;
        }

        return 1;
    }

    private int add(byte obsBase1, byte qual1, byte obsBase2, byte qual2) {
        return add(obsBase1, qual1, obsBase2, qual2, 1);
    }

    // -------------------------------------------------------------------------------------
    //
    // Dealing with the cache routines
    //
    // -------------------------------------------------------------------------------------

    static DiploidSNPGenotypeLikelihoods[][][][][] CACHE = new DiploidSNPGenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_QUAL_SCORE+1][BaseUtils.BASES.length+1][QualityUtils.MAX_QUAL_SCORE+1][MAX_PLOIDY];

    protected boolean inCache(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        return getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy) != null;
    }

    protected DiploidSNPGenotypeLikelihoods getCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        DiploidSNPGenotypeLikelihoods gl = getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy);
        if ( gl == null )
            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base1=%c, qual1=%d, base2=%c, qual2=%d, ploidy=%d",
                    observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy));
        return gl;
    }

    protected DiploidSNPGenotypeLikelihoods calculateCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        DiploidSNPGenotypeLikelihoods gl = calculateGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);
        setCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy, gl);
        return gl;
    }

    protected void setCache( DiploidSNPGenotypeLikelihoods[][][][][] cache,
                             byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy,
                             DiploidSNPGenotypeLikelihoods val ) {
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        int j = qualityScore1;
        int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        int l = qualityScore2;
        int m = ploidy;

        cache[i][j][k][l][m] = val;
    }

    protected DiploidSNPGenotypeLikelihoods getCache(DiploidSNPGenotypeLikelihoods[][][][][] cache,
                                            byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
        int j = qualityScore1;
        int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
        int l = qualityScore2;
        int m = ploidy;
        return cache[i][j][k][l][m];
    }

    protected DiploidSNPGenotypeLikelihoods calculateGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
        double[] log10FourBaseLikelihoods = computeLog10Likelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);

        try {

            DiploidSNPGenotypeLikelihoods gl = (DiploidSNPGenotypeLikelihoods)this.clone();
            gl.setToZero();

            // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
            for ( DiploidGenotype g : DiploidGenotype.values() ) {

                // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
                double p_base = 0.0;
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);
                double likelihood = log10(p_base);

                gl.log10Likelihoods[g.ordinal()] += likelihood;
                gl.log10Posteriors[g.ordinal()] += likelihood;
            }

            if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
            }

            return gl;

         } catch ( CloneNotSupportedException e ) {
             throw new RuntimeException(e);
         }
    }

    /**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase1  the base observed on the 1st read of the fragment
     * @param qualityScore1  the qual of the base on the 1st read of the fragment, or zero if NA
     * @param observedBase2  the base observed on the 2nd read of the fragment
     * @param qualityScore2  the qual of the base on the 2nd read of the fragment, or zero if NA
     * @return likelihoods for this observation or null if the base was not considered good enough to add to the likelihoods (Q0 or 'N', for example)
     */
    protected double[] computeLog10Likelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
        double[] log10FourBaseLikelihoods = baseZeros.clone();

        for ( byte trueBase : BaseUtils.BASES ) {
            double likelihood = 0.0;

            for ( byte fragmentBase : BaseUtils.BASES ) {
                double log10FragmentLikelihood = (trueBase == fragmentBase ? log10_1_minus_PCR_error : log10_PCR_error_3);
                if ( qualityScore1 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase1, fragmentBase, qualityScore1);
                }
                if ( qualityScore2 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase2, fragmentBase, qualityScore2);
                }

                //if ( VERBOSE ) {
                //    System.out.printf("  L(%c | b=%s, Q=%d) = %f / %f%n",
                //            observedBase, trueBase, qualityScore, pow(10,likelihood) * 100, likelihood);
                //}

                likelihood += pow(10, log10FragmentLikelihood);
            }

            log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(trueBase)] = log10(likelihood);
        }

        return log10FourBaseLikelihoods;
    }

    /**
     *
     * @param observedBase observed base
     * @param chromBase    target base
     * @param qual         base quality
     * @return log10 likelihood
     */
    protected double log10PofObservingBaseGivenChromosome(byte observedBase, byte chromBase, byte qual) {

        double logP;

        if ( observedBase == chromBase ) {
            // the base is consistent with the chromosome -- it's 1 - e
            //logP = oneMinusData[qual];
            double e = pow(10, (qual / -10.0));
            logP = log10(1.0 - e);
        } else {
            // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
            logP = qual / -10.0 + (-log10_3);
        }

        //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
        return logP;
    }

    /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p
     * @param ignoreBadBases
     * @param capBaseQualsAtMappingQual
     * @param minBaseQual
     * @return
     */
    private static byte qualToUse(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase( p.getBase() ) ) {
            return 0;
        } else {
            byte qual = p.getQual();

            if ( qual > SAMUtils.MAX_PHRED_SCORE )
                throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
            if ( capBaseQualsAtMappingQual )
                qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
            if ( (int)qual < minBaseQual )
                qual = (byte)0;

            return qual;
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
     * Return a string representation of this object in a moderately usable form
     *
     * @return string representation
     */
    public String toString() {
        double sum = 0;
        StringBuilder s = new StringBuilder();
        for (DiploidGenotype g : DiploidGenotype.values()) {
            s.append(String.format("%s %.10f ", g, log10Likelihoods[g.ordinal()]));
			sum += Math.pow(10,log10Likelihoods[g.ordinal()]);
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
                if ( ! MathUtils.wellFormedDouble(log10Likelihoods[i]) || ! MathUtils.isNegativeOrZero(log10Likelihoods[i]) ) {
                    bad = String.format("Likelihood %f is badly formed", log10Likelihoods[i]);
                } else if ( ! MathUtils.wellFormedDouble(log10Posteriors[i]) || ! MathUtils.isNegativeOrZero(log10Posteriors[i]) ) {
                    bad = String.format("Posterior %f is badly formed", log10Posteriors[i]);
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

    //
    // Constant static data
    //
    private final static double[] genotypeZeros = new double[DiploidGenotype.values().length];
    private final static double[] baseZeros = new double[BaseUtils.BASES.length];

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            genotypeZeros[g.ordinal()] = 0.0;
        }
        for ( byte base : BaseUtils.BASES ) {
            baseZeros[BaseUtils.simpleBaseToBaseIndex(base)] = 0.0;
        }
    }
}
