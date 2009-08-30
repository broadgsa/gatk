package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.BaseTest;
import org.junit.Ignore;

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
 * for homozygous genotypes and
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
@Ignore
public class NewHotnessGenotypeLikelihoodsTest extends BaseTest {
    int x;
/*    private int coverage = 0;
    private double[] likelihoods = null;
    private double[] priors = null;
    private double[] posteriors = null;

    GenotypeLikelihoods();
    GenotypeLikelihoods(char ref, double heterozygosity)
    GenotypeLikelihoods(char ref, double priorHomRef, double priorHet, double priorHomVar);
    GenotypeLikelihoods(double[] log10Priors)
    double[] getLikelihoods()
    double getLikelihood(DiploidGenotype g);
    double[] getPosteriors();
    double getPosterior(DiploidGenotype g);

    double[] getPriors()
    double getPrior(DiploidGenotype g)
    int getCoverage()
    boolean isFilteringQ0Bases()
    void filterQ0Bases(boolean filterQ0Bases)
    int add(char ref, char read, byte qual)
    int add(char observedBase, byte qualityScore)
    private boolean badBase(char observedBase)
    int add(ReadBackedPileup pileup, boolean ignoreBadBases)
    private double calculateBaseLikelihood(char read, String genotype, byte qual)
    String toString()
    boolean validate(boolean throwException)*/
}