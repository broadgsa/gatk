package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;

public class IVFSecondaryBases implements IndependentVariantFeature {
    private double[] p2on = { 0.000, 0.302, 0.366, 0.142, 0.000, 0.548, 0.370, 0.000, 0.319, 0.000 };
    private double[] p2off = { 0.480, 0.769, 0.744, 0.538, 0.575, 0.727, 0.768, 0.589, 0.762, 0.505 };

    /**
     * Method so that features can initialize themselves based on a short argument string. At the moment, each feature is
     * responsible for interpreting their own argument string.
     *
     * @param arguments the arguments!
     */
    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(";");
            String[] argOnPieces = argPieces[0].split(",");
            String[] argOffPieces = argPieces[0].split(",");

            for (int genotypeIndex = 0; genotypeIndex < 10; genotypeIndex++) {
                p2on[genotypeIndex] = Double.valueOf(argOnPieces[genotypeIndex]);
                p2off[genotypeIndex] = Double.valueOf(argOffPieces[genotypeIndex]);
            }
        }
    }

    /**
     * Method to compute the result of this feature for each of the ten genotypes.  The return value must be a double array
     * of length 10 (one for each genotype) and the value must be in log10-space.
     *
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return a ten-element array of log-likelihood result of the feature applied to each genotype
     */
    public double[] compute(char ref, LocusContext context) {
        double[] likelihoods = new double[10];

        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String primaryBases = pileup.getBases();
        String secondaryBases = pileup.getSecondaryBasePileup();

        for (int genotypeIndex = 0; genotypeIndex < Genotype.values().length; genotypeIndex++) {
            char firstAllele = Genotype.values()[genotypeIndex].toString().charAt(0);
            char secondAllele = Genotype.values()[genotypeIndex].toString().charAt(1);

            int offIsGenotypic = 0;
            int offTotal = 0;

            int onIsGenotypic = 0;
            int onTotal = 0;

            for (int pileupIndex = 0; pileupIndex < primaryBases.length(); pileupIndex++) {
                char primaryBase = primaryBases.charAt(pileupIndex);

                if (secondaryBases != null) {
                    char secondaryBase = secondaryBases.charAt(pileupIndex);

                    if (primaryBase != firstAllele && primaryBase != secondAllele) {
                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
                            offIsGenotypic++;
                        }
                        offTotal++;
                    } else {
                        if (secondaryBase == firstAllele || secondaryBase == secondAllele) {
                            onIsGenotypic++;
                        }
                        onTotal++;
                    }
                }
            }

            double offPrior = MathUtils.binomialProbability(offIsGenotypic, offTotal, p2off[genotypeIndex]);
            double onPrior = MathUtils.binomialProbability(onIsGenotypic, onTotal, p2on[genotypeIndex]);

            double logOffPrior = MathUtils.compareDoubles(offPrior, 0.0, 1e-10) == 0 ? Math.log10(Double.MIN_VALUE) : Math.log10(offPrior);
            double logOnPrior = MathUtils.compareDoubles(onPrior, 0.0, 1e-10) == 0 ? Math.log10(Double.MIN_VALUE) : Math.log10(onPrior);

            likelihoods[genotypeIndex] += logOffPrior + logOnPrior;
        }

        return likelihoods;
    }
}
