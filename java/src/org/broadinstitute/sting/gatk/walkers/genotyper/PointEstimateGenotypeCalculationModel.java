package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import java.util.*;

public class PointEstimateGenotypeCalculationModel extends EMGenotypeCalculationModel {

    protected PointEstimateGenotypeCalculationModel() {}

    // overload this method so we can special-case the single sample
    public List<GenotypeCall> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // we don't actually want to run EM for single samples
        if ( samples.size() == 1 ) {

            // split the context (so we can get forward/reverse contexts for free)
            HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context, new int[4]);
            if ( contexts == null )
                return null;

            // get the context for the sample
            String sample = samples.iterator().next();
            AlignmentContextBySample sampleContext = contexts.get(sample);

            // if there were no good bases, the context wouldn't exist
            if ( sampleContext == null )
                return null;

            callsMetrics.nCalledBases++;

            // create the genotype call object
            Pair<ReadBackedPileup, GenotypeLikelihoods> discoveryGL = getSingleSampleLikelihoods(ref, sampleContext, priors, StratifiedContext.OVERALL);
            GenotypeCall call = new GenotypeCall(sample, context.getLocation(), ref, discoveryGL.second, discoveryGL.first);

            if ( GENOTYPE_MODE || call.isVariant(call.getReference()) ) {
                double confidence = (GENOTYPE_MODE ? call.getNegLog10PError() : call.toVariation().getNegLog10PError());
                if ( confidence >= LOD_THRESHOLD ) {
                    callsMetrics.nConfidentCalls++;
                    out.addGenotypeCall(call);
                } else {
                    callsMetrics.nNonConfidentCalls++;
                }
            }
            return Arrays.asList(call);
        }

        return super.calculateGenotype(tracker, ref, context, priors);
    }

    private Pair<ReadBackedPileup, GenotypeLikelihoods> getSingleSampleLikelihoods(char ref, AlignmentContextBySample sampleContext, DiploidGenotypePriors priors, StratifiedContext contextType) {
        // create the pileup
        AlignmentContext myContext = sampleContext.getContext(contextType);
        ReadBackedPileup pileup = new ReadBackedPileup(ref, myContext);
        pileup.setIncludeDeletionsInPileupString(true);

        // create the GenotypeLikelihoods object
        GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
        GL.setVerbose(VERBOSE);
        GL.add(pileup, true);
        return new Pair<ReadBackedPileup, GenotypeLikelihoods>(pileup, GL);
    }

    protected double[] initializeAlleleFrequencies(int numSamplesInContext,  int[] baseCounts) {
        // An intelligent guess would be the observed base ratios
        // (plus some small number to account for sampling issues).
        int totalObservedBases = 0;
        for (int i = 0; i < 4; i++)
            totalObservedBases += baseCounts[i];

        double[] alleleFrequencies = new double[4];
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = ((double)baseCounts[i] + DiploidGenotypePriors.HUMAN_HETEROZYGOSITY) / (double)totalObservedBases;
        return alleleFrequencies;
    }

    protected HashMap<String, GenotypeLikelihoods> initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, double[] alleleFrequencies, DiploidGenotypePriors priors, StratifiedContext contextType) {
        HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();

        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, AFPriors, defaultPlatform);
            GL.setVerbose(VERBOSE);
            GL.add(pileup, true);

            GLs.put(sample, GL);
        }

        return GLs;
    }

    protected double[] calculateAlleleFrequencyPosteriors(HashMap<String, GenotypeLikelihoods> GLs) {
        double[] newAlleleLikelihoods = new double[4];

        for ( GenotypeLikelihoods GL : GLs.values() ) {
            double[] normalizedPosteriors = GL.getNormalizedPosteriors();

            // calculate the posterior weighted frequencies for this sample
            double[] personalAllelePosteriors = new double[4];
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                double posterior = normalizedPosteriors[g.ordinal()] / 2.0;   // each base gets half the probability
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base1)] += posterior;
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base2)] += posterior;
            }

            for (int i = 0; i < 4; i++)
                newAlleleLikelihoods[i] += personalAllelePosteriors[i];
        }

        // normalize
        double sum = 0;
        for (int i = 0; i < 4; i++)
            sum += newAlleleLikelihoods[i];
        for (int i = 0; i < 4; i++)
            newAlleleLikelihoods[i] /= sum;

        return newAlleleLikelihoods;
    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods(HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies) {
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);
        for ( GenotypeLikelihoods GL : GLs.values() )
            GL.setPriors(AFPriors);
    }

    protected EMOutput computePofF(char ref, HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies, int numSamplesInContext) {

        // compute pD and pNull without allele frequencies
        double pD = compute_pD(GLs);
        double pNull = compute_pNull(ref, GLs);
        logger.debug("Original pD=" + pD + ", pNull=" + pNull);

        // compute p0
        double pVar = 0.0;
        for (int i = 1; i < numSamplesInContext; i++)
            pVar += heterozygosity/(double)i;
        double p0 = Math.log10(1.0 - pVar);

        // compute actual priors: theta / MAF
        double MAF;
        Integer[] sortedIndexes = Utils.SortPermutation(alleleFrequencies);
        if ( sortedIndexes[3] != BaseUtils.simpleBaseToBaseIndex(ref) )
            MAF = alleleFrequencies[sortedIndexes[3]];
        else
            MAF = alleleFrequencies[sortedIndexes[2]];

        //  compute pF
        double pF;
        double expectedChromosomes = 2.0 * (double)numSamplesInContext * MAF;
        if ( expectedChromosomes < 1.0 )
            pF = p0;
        else
            pF = Math.log10(heterozygosity / expectedChromosomes);
        logger.debug("p0=" + p0 + ", pF=" + pF);

        pD += pF;
        pNull += p0;
        logger.debug("Final pD=" + pD + ", pNull=" + pNull);

        return new EMOutput(pD, pNull, pF, GLs);
    }

    private double compute_pD(HashMap<String, GenotypeLikelihoods> GLs) {
        double pD = 0.0;
        for ( GenotypeLikelihoods GL : GLs.values() ) {
            double sum = 0.0;
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                sum += Math.pow(10, GL.getPosterior(g));
            }
            pD += Math.log10(sum);
        }
        return pD;
    }

    private double compute_pNull(char ref, HashMap<String, GenotypeLikelihoods> GLs) {
        // compute null likelihoods
        double[] alleleLikelihoods = new double[4];
        for (int i = 0; i < 4; i++)
            alleleLikelihoods[i] = 1e-6/3.0;
        alleleLikelihoods[BaseUtils.simpleBaseToBaseIndex(ref)] = 1.0-1e-6;
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleLikelihoods);

        HashMap<String, GenotypeLikelihoods> GL_null = new HashMap<String, GenotypeLikelihoods>();
        try {
            for ( String sample : GLs.keySet() ) {
                GenotypeLikelihoods GL = (GenotypeLikelihoods)GLs.get(sample).clone();
                GL.setPriors(AFPriors);
                GL_null.put(sample, GL);
            }
        } catch (CloneNotSupportedException e) {
            throw new StingException("Clone() not supported for given GenotypeLikelihoods subclass?");
        }

        return compute_pD(GL_null);
    }
}