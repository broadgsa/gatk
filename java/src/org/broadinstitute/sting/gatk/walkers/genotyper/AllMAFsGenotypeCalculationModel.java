package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public class AllMAFsGenotypeCalculationModel extends EMGenotypeCalculationModel {

    protected AllMAFsGenotypeCalculationModel() {}

    // because the null allele frequencies are constant for a given N,
    // we cache the results to avoid having to recompute everything
    private HashMap<Integer, double[]> nullAlleleFrequencyCache = new HashMap<Integer, double[]>();

    // the allele frequencies
    private double[][] alleleFrequencies = new double[3][];
    private double[][] oldAlleleFrequencies;

    // keep track of whether or not a given MAF is stable
    private boolean[] frequencyStabilityArray = new boolean[3];

    // the minimum and actual number of points in our allele frequency estimation
    private static final int MIN_ESTIMATION_POINTS = 1000;
    private int estimationPoints;

    // the GenotypeLikelihoods map
    private HashMap<String, AlleleSpecificGenotypeLikelihoods> GLs = new HashMap<String, AlleleSpecificGenotypeLikelihoods>();


    protected void initializeAlleleFrequencies(int numSamples, char ref) {
        // first, initialize the stability array to "unstable"
        for (int i = 0; i < 3; i++)
            frequencyStabilityArray[i] = false;

        // calculate the number of estimation points to use:
        // it's either MIN_ESTIMATION_POINTS or 2N if that's larger
        // (add 1 for allele frequency of zero)
        estimationPoints = Math.max(MIN_ESTIMATION_POINTS, 2 * numSamples) + 1;

        for (int alt = 0; alt < 3; alt++)
            alleleFrequencies[alt] = getNullAlleleFrequencies(estimationPoints);

        for (int i = 1; i < estimationPoints; i++)
            logger.debug("Initial allele frequency for MAF=" + i + ": " + alleleFrequencies[0][i]);
    }

    private double[] getNullAlleleFrequencies(int N) {
        double[] AFs = nullAlleleFrequencyCache.get(N);

        // if it hasn't been calculated yet, do so now
        if ( AFs == null ) {

            // calculate sum(1/i)
            double denominator = 0.0;
            for (int i = 1; i < N; i++)
                denominator += 1.0 / (double)i;

            // set up delta
            double delta = 1.0 / denominator;

            // calculate the null allele frequencies
            AFs = new double[N];
            for (int i = 1; i < N; i++)
                AFs[i] = Math.log10(delta / (double)i);

            nullAlleleFrequencyCache.put(N, AFs);
        }

        return AFs.clone();
    }

    protected void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {
        GLs.clear();

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.setVerbose(VERBOSE);
            GL.add(pileup, true);

            GLs.put(sample, new AlleleSpecificGenotypeLikelihoods(ref, GL));
        }
    }

    protected void calculateAlleleFrequencyPosteriors() {

    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods() {

    }

    protected boolean isStable() {
        // We consider the EM stable when the MAF doesn't change more than EM_STABILITY_METRIC
        // We compute this separately for all of the alternate alleles
        for (int i = 0; i < 3; i++) {
            // if we've already determined that a MAF is stable, don't recalculate
            if ( frequencyStabilityArray[i] )
                continue;

            // determine change
            double AF_delta = 0.0;
            for (int j = 1; j < estimationPoints; j++)
                AF_delta += Math.abs(oldAlleleFrequencies[i][j] - alleleFrequencies[i][j]);

            // if it's not stable, we're done
            if (AF_delta > EM_STABILITY_METRIC)
                return false;

            // it's stable, so record that fact in the stability array
            frequencyStabilityArray[i] = true;
        }

        // if we got here, then we're stable
        return true;
    }

    protected EMOutput computePofF(char ref) {
        return null;
    }


    /**
     * A class for the allele-specific genotype likelihoods
    */
    protected class AlleleSpecificGenotypeLikelihoods {
        private HashMap<String, Pair<double[], double[]>> GLs = new HashMap<String, Pair<double[], double[]>>();

        AlleleSpecificGenotypeLikelihoods(char ref, GenotypeLikelihoods GL) {
            double[] likelihoods = GL.getLikelihoods();
            double[] posteriors = GL.getPosteriors();



            // get the ref likelihood/posterior
            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
            double refLikelihood = GL.getLikelihood(refGenotype);
            double refPosterior = GL.getPosterior(refGenotype);
            String refStr = String.valueOf(ref);

            for ( char base : BaseUtils.BASES ) {
                if ( base == ref )
                    continue;

                // get hom var and het likelihoods
                double homLikelihood = GL.getLikelihood(DiploidGenotype.createHomGenotype(base));
                double hetLikelihood = GL.getLikelihood(DiploidGenotype.valueOf(refStr + String.valueOf(base)));


            }
        }




    }
}