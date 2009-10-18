package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;

import java.util.*;

public class AllMAFsGenotypeCalculationModel extends EMGenotypeCalculationModel {

    protected AllMAFsGenotypeCalculationModel() {}

    private double[] alleleFrequencies;



    protected void initializeAlleleFrequencies(int numSamples, char ref) {
        // we have 2N possible allele frequencies in pileup
        int possibleMAFs = 2 * numSamples;
        alleleFrequencies = new double[possibleMAFs];

        // calculate sum(1/i) for i from 1 to 2N
        double denominator = 0.0;
        for (int i = 1; i <= possibleMAFs; i++)
            denominator += 1.0 / (double)i;

        // set up delta
        double delta = 1.0 / denominator;

        // calculate the null allele frequencies
        for (int i = 1; i <= possibleMAFs; i++)
            alleleFrequencies[i-1] = Math.log10(delta / (double)i);

        for (int i = 0; i < possibleMAFs; i++)
            logger.debug("Initial allele frequency for MAF=" + (i+1) + ": " + alleleFrequencies[i]);
    }

    protected void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {
        HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.setVerbose(VERBOSE);
            GL.add(pileup, true);

            GLs.put(sample, GL);
        }
    }

    protected void calculateAlleleFrequencyPosteriors() {

    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods() {

    }

    protected boolean isStable() {
        return true;
    }

    protected EMOutput computePofF(char ref) {
        return null;
    }
}