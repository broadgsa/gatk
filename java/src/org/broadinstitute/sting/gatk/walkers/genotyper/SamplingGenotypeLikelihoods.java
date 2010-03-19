package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.util.Arrays;
import java.util.Comparator;

public class SamplingGenotypeLikelihoods extends GenotypeLikelihoods {
    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     * @param m base model
     */
    public SamplingGenotypeLikelihoods(BaseMismatchModel m) {
        super(m);
        enableCacheFlag = false;
    }

    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     * @param m base model
     * @param pl default platform
     */
    public SamplingGenotypeLikelihoods(BaseMismatchModel m, EmpiricalSubstitutionProbabilities.SequencerPlatform pl) {
        super(m, pl);
        enableCacheFlag = false;
    }

    /**
     * Create a new GenotypeLikelhoods object with given priors for each diploid genotype
     *
     * @param m base model
     * @param priors priors
     */
    public SamplingGenotypeLikelihoods(BaseMismatchModel m, DiploidGenotypePriors priors) {
        super(m, priors);
        enableCacheFlag = false;
    }

    /**
     * Create a new GenotypeLikelhoods object with given priors for each diploid genotype
     *
     * @param m base model
     * @param priors priors
     * @param pl default platform
     */
    public SamplingGenotypeLikelihoods(BaseMismatchModel m, DiploidGenotypePriors priors, EmpiricalSubstitutionProbabilities.SequencerPlatform pl) {
        super(m, priors, pl);
        enableCacheFlag = false;
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    /**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup         read pileup
     * @param ignoreBadBases should we ignore bad bases
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases) {
        // we're actually going to be happy with our homozygous theories
        int n = super.add(pileup, ignoreBadBases);

        // Now, loop over the heterygous theories and do our fancy sampling calculation
        FourBaseProbabilities[] probs = fourBaseProbVector(n, pileup, ignoreBadBases);
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            if ( g.isHet() ) {
                //System.out.printf("DEBUG: computing best k for %s at %s%n", g, pileup.getLocation());
                double likelihood = calculateSamplingHetGenotypeLikelihood(probs, g, n);
                //System.out.printf("DEBUG: L was %.2f%n", likelihood);
                if ( likelihood != 1 ) {
                    log10Likelihoods[g.ordinal()] = likelihood;
                    log10Posteriors[g.ordinal()] = likelihood + priors.getPrior(g);
                }
            }
        }

        return n;
    }

    public int add(char observedBase, byte qualityScore, SAMRecord read, int offset) {
        throw new UnsupportedOperationException("BUG: Sampling genotype likelihoods does not support sequential addition of bases; use add(ReadBackedPileup) instead");
    }


    private double calculateSamplingHetGenotypeLikelihood(FourBaseProbabilities[] probs, DiploidGenotype g, int goodBaseDepth) {
        double bestLog10L = 1;
        //int bestK = -1;

        for ( int k = 0; k < goodBaseDepth; k++ ) {
            // for every possible sampling depth of the A allele
            double log10BinomialSampling = Math.log10(MathUtils.binomialProbabilityLog(k, goodBaseDepth, 0.5));
            double log10Het = mostLikelyKBases(probs, k, BaseUtils.simpleBaseToBaseIndex(g.base1), BaseUtils.simpleBaseToBaseIndex(g.base2));

            double log10L = log10BinomialSampling + log10Het;
            //System.out.printf("  DEBUG: computing best k = %d, bin = %.2f, het = %.2f, log10L = %.2f%n", k, log10BinomialSampling, log10Het, log10L);

            if ( log10L > bestLog10L || bestLog10L == 1 ) {
                bestLog10L = log10L;
                //bestK = k;
            }
        }

        //System.out.printf("  DEBUG: best k = %d with log10L of %.2f%n", bestK, bestLog10L);
        return bestLog10L;
    }

    FourBaseProbabilities[] fourBaseProbVector(int n, ReadBackedPileup pileup, boolean ignoreBadBases) {
        FourBaseProbabilities[] probs = new FourBaseProbabilities[n];

        int i = 0;
        for ( PileupElement p : pileup ) {
            char observedBase = (char)p.getBase();
            byte qualityScore = p.getQual();
            SAMRecord read = p.getRead();
            int offset = p.getOffset();

            if ( ! ignoreBadBases || ! badBase(observedBase) ) {
                FourBaseProbabilities fbl = fourBaseLikelihoods.computeLog10Likelihoods(observedBase, qualityScore, read, offset);

                if ( fbl != null ) {
                    probs[i++] = fbl;
                }
            }
        }

        return probs;
    }

    class FourBaseComparator implements Comparator<FourBaseProbabilities> {
        int index = 0;
        public FourBaseComparator(int i) { index = i; }
        public int compare(FourBaseProbabilities a, FourBaseProbabilities b) {
            return -1 * Double.compare(a.getLog10Likelihood(index), b.getLog10Likelihood(index));
        }
    }

    protected static final double log103 = log10(3.0);
    double mostLikelyKBases(FourBaseProbabilities[] probs, int k, int base1, int base2) {
        if ( k == 0 ) // optimization -- if k > 1, we've already sorted the vector
            Arrays.sort(probs, new FourBaseComparator(base1));

        double log10L = 0;
        for ( int i = 0; i < probs.length; i++ ) {
            int base = i < k ? base1 : base2;
            log10L += probs[i].getLog10Likelihood(base) - log103;
            //System.out.printf("%3d %d %.2f %.2f %.2f%n", i, base, probs[i].getLog10Likelihood(base), probs[i].getLog10Likelihood(base1), probs[i].getLog10Likelihood(base2));
        }

        return log10L;
    }
}