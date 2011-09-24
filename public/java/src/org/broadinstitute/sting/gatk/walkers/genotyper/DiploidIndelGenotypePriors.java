package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.walkers.indels.HaplotypeIndelErrorModel;
import org.broadinstitute.sting.utils.MathUtils;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: Sep 30, 2010
 * Time: 1:47:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class DiploidIndelGenotypePriors implements GenotypePriors {
    // --------------------------------------------------------------------------------------------------------------
    //
    // Constants and static information
    //
    // --------------------------------------------------------------------------------------------------------------
    public static final double INDEL_HETEROZYGOSITY = 1e-4;

    private final static double[] flatPriors = new double[DiploidGenotype.values().length];

    // --------------------------------------------------------------------------------------------------------------
    //
    // Diploid priors
    //
    // --------------------------------------------------------------------------------------------------------------
    private double[] priors = null;

    /**
     * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
     */
    public DiploidIndelGenotypePriors() {
        priors = flatPriors.clone();
    }

    public DiploidIndelGenotypePriors(double indelHeterozygosity, int eventLength, int haplotypeSize)  {
        double varPrior = getHaplotypePriors(indelHeterozygosity, eventLength, haplotypeSize);
        priors[2] = Math.log10(varPrior*varPrior);
        priors[1] = Math.log10(2*varPrior*(1-varPrior));
        priors[0] = Math.log10((1-varPrior)*(1-varPrior));
 
    }


     /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors;
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    public double getHeterozygosity() { return INDEL_HETEROZYGOSITY; }

    public boolean validate(boolean throwException) {
        try {

            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(priors[i]) || ! MathUtils.isNegativeOrZero(priors[i]) ) {
                    String bad = String.format("Prior %f is badly formed %b", priors[i], MathUtils.isNegativeOrZero(priors[i]));
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

    public double getHaplotypePriors(double indelHeterozygosity, int eventLength, int haplotypeSize) {
        // compute prior likelihoods on haplotypes.
        // In general, we'll assume: even spread of indels throughout genome (not true, but simplifying assumption),
        // and memoryless spread (i.e. probability that an indel lies in an interval A is independent of probability of
        // another indel lying in interval B iff A and B don't overlap), then we can approximate inter-indel distances
        // by an exponential distribution of mean 1/theta (theta = heterozygozity), and the number of indels on an interval
        // of size L is Poisson-distributed with parameter lambda = theta*L.

        // Since typically, for small haplotype sizes and human heterozygozity, lambda will be <<1, we'll further approximate it
        // by assuming that only one indel can happen in a particular interval, with Pr(indel present) = lambda*exp(-lambda), and
        // pr(no indel) = 1-lambda*exp(-lambda) ~= exp(-lambda) for small lambda.

        // We also assume that a deletion is equally likely as an insertion (empirical observation, see e.g. Mills et al, Genome Research 2006)
        // and we assume the following frequency spectrum for indel sizes Pr(event Length = L)= K*abs(L)^(-1.89)*10^(-0.015*abs(L)),
        // taking positive L = insertions, negative L = deletions. K turns out to be about 1.5716 for probabilities to sum to one.
        // so -10*log10(Pr event Length = L) =-10*log10(K)+ 18.9*log10(abs(L)) + 0.15*abs(L).
        // Hence, Pr(observe event size = L in interval) ~ Pr(observe event L | event present) Pr (event present in interval)
        // and -10*log10(above) = -10*log10(K)+ 18.9*log10(abs(L)) + 0.15*abs(L) - 10*log10(theta*L), and we ignore terms that would be
        // added to ref hypothesis.
        // Equation above is prior model.

        double lambda = (double)haplotypeSize * indelHeterozygosity;
        return HaplotypeIndelErrorModel.probToQual(lambda)-HaplotypeIndelErrorModel.probToQual(eventLength)*1.89 + 0.15*eventLength
                + HaplotypeIndelErrorModel.probToQual(1.5716)+ HaplotypeIndelErrorModel.probToQual(0.5);



    }


    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
        }
    }
}

