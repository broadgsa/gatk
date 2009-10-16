package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;

import java.util.*;

public class AllMAFsGenotypeCalculationModel extends EMGenotypeCalculationModel {


    protected AllMAFsGenotypeCalculationModel() {}

    protected double[] initializeAlleleFrequencies(int numSamples, char ref) {
        double[] alleleFrequencies = new double[2 * numSamples + 1];   // 2N + 1 possible minor alleles



        return alleleFrequencies;
    }

    protected HashMap<String, GenotypeLikelihoods> initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, double[] alleleFrequencies, DiploidGenotypePriors priors, StratifiedContext contextType) {
        return null;
    }

    protected double[] calculateAlleleFrequencyPosteriors(HashMap<String, GenotypeLikelihoods> GLs) {
        return null;
    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods(HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies) {

    }

    protected EMOutput computePofF(char ref, HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies, int numSamplesInContext) {
        return null;
    }
}