package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public class PointEstimateGenotypeCalculationModel extends EMGenotypeCalculationModel {

    protected PointEstimateGenotypeCalculationModel() {}

    // the allele frequencies
    private double[] alleleFrequencies = new double[4];
    private double[] oldAlleleFrequencies;

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();

    // Allele frequency initialization values from the original MSG code (so we can be consistent)
    private static final double NON_REF = 0.0005002502;  // heterozygosity / (2 * sqrt(1-heterozygosity)
    private static final double REF = 0.9994999;         //sqrt(1-heterozygosity)


    // overload this method so we can special-case the single sample
    public Pair<List<Genotype>, GenotypeMetaData> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // we don't actually want to run EM for single samples
        if ( samples.size() == 1 ) {

            // split the context (so we can get forward/reverse contexts for free)
            HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context);
            if ( contexts == null )
                return null;

            // get the context for the sample
            String sample = samples.iterator().next();
            AlignmentContextBySample sampleContext = contexts.get(sample);

            // if there were no good bases, the context wouldn't exist
            if ( sampleContext == null )
                return null;

            // create the genotype call object
            // create the call
            Genotype call = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT);
            call.setReference(ref);
            call.setLocation(context.getLocation());

            Pair<ReadBackedPileup, GenotypeLikelihoods> discoveryGL = getSingleSampleLikelihoods(ref, sampleContext, priors, StratifiedContext.OVERALL);

            if ( call instanceof ReadBacked ) {
                ((ReadBacked)call).setReads(discoveryGL.first.getReads());
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(discoveryGL.second.getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(discoveryGL.second.getPosteriors());
            }

            return new Pair<List<Genotype>, GenotypeMetaData>(Arrays.asList(call), null);
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
        GL.add(pileup, true);
        return new Pair<ReadBackedPileup, GenotypeLikelihoods>(pileup, GL);
    }

    protected void initializeAlleleFrequencies(int numSamplesInContext, char ref) {
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = NON_REF;
        alleleFrequencies[BaseUtils.simpleBaseToBaseIndex(ref)] = REF;

        for (int i = 0; i < 4; i++)
            logger.debug("Initial allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + alleleFrequencies[i]);
    }

    protected void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {
        GLs.clear();

        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, AFPriors, defaultPlatform);
            GL.add(pileup, true);

            GLs.put(sample, GL);
        }
    }

    private static DiploidGenotypePriors calculateAlleleFreqBasedPriors(double[] alleleFreqs) {
        // convert to log-space
        double[] log10Freqs = new double[4];
        for (int i = 0; i < 4; i++)
            log10Freqs[i] = Math.log10(alleleFreqs[i]);

        double[] alleleFreqPriors = new double[10];

        // this is the Hardy-Weinberg based allele frequency (p^2, q^2, 2pq)
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            alleleFreqPriors[g.ordinal()] = log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base1)] + log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base2)];
            // add a factor of 2 for the 2pq case
            if ( g.isHet() )
                alleleFreqPriors[g.ordinal()] += Math.log10(2);
        }

        return new DiploidGenotypePriors(alleleFreqPriors);
    }

    protected void calculateAlleleFrequencyPosteriors() {
        // initialization
        oldAlleleFrequencies = alleleFrequencies.clone();
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = 0.0;

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
                alleleFrequencies[i] += personalAllelePosteriors[i];
        }

        // normalize
        double sum = 0.0;
        for (int i = 0; i < 4; i++)
            sum += alleleFrequencies[i];
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] /= sum;

        for (int i = 0; i < 4; i++)
            logger.debug("New allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + alleleFrequencies[i]);
    }

    protected void applyAlleleFrequencyToGenotypeLikelihoods() {
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);
        for ( GenotypeLikelihoods GL : GLs.values() )
            GL.setPriors(AFPriors);
    }

    protected boolean isStable() {
        // We consider the EM stable when the MAF doesn't change more than EM_STABILITY_METRIC
        double AF_delta = 0.0;
        for (int i = 0; i < 4; i++)
            AF_delta += Math.abs(oldAlleleFrequencies[i] - alleleFrequencies[i]);

        return (AF_delta < EM_STABILITY_METRIC);
    }

    protected EMOutput computePofF(char ref) {
        // some debugging output
        for ( String sample : GLs.keySet() )
            logger.debug("GenotypeLikelihoods for sample " + sample + ": " + GLs.get(sample).toString());

        // compute pD and pNull without allele frequencies
        double pD = compute_pD(GLs);
        double pNull = compute_pNull(ref, GLs);
        logger.debug("Original pD=" + pD + ", pNull=" + pNull);

        // compute p0
        double pVar = 0.0;
        for (int i = 1; i < GLs.size(); i++)
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
        double expectedChromosomes = 2.0 * (double)GLs.size() * MAF;
        if ( expectedChromosomes < 1.0 )
            pF = p0;
        else
            pF = Math.log10(heterozygosity / expectedChromosomes);
        logger.debug("p0=" + p0 + ", pF=" + pF);

        pD += pF;
        pNull += p0;
        logger.debug("Final pD=" + pD + ", pNull=" + pNull);

        return new EMOutput(pD, pNull, pF, MAF, GLs);
    }

    private static double compute_pD(HashMap<String, GenotypeLikelihoods> GLs) {
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

    private static double compute_pNull(char ref, HashMap<String, GenotypeLikelihoods> GLs) {
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