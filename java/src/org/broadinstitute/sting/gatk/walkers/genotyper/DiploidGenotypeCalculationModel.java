package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

public class DiploidGenotypeCalculationModel extends JointEstimateGenotypeCalculationModel {

    protected DiploidGenotypeCalculationModel() {}

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();


    protected HashMap<String, AlignmentContextBySample> createContexts(AlignmentContext context) {
        return splitContextBySample(context);
    }

    protected int getNSamples(HashMap<String, AlignmentContextBySample> contexts) {
        return contexts.size();
    }

    protected void initializeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        GLs.clear();

        // use flat priors for GLs
        DiploidGenotypePriors priors = new DiploidGenotypePriors();

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = new GenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.add(pileup, true);
            GLs.put(sample, GL);
        }
    }

    protected double computeLog10PofDgivenAFi(DiploidGenotype refGenotype, DiploidGenotype hetGenotype, DiploidGenotype homGenotype, double f) {
        double PofDgivenAFi = 0.0;

        // for each sample
        for ( GenotypeLikelihoods GL : GLs.values() ) {

            double[] posteriors = GL.getPosteriors();

            double[] allelePosteriors = new double[] { posteriors[refGenotype.ordinal()], posteriors[hetGenotype.ordinal()], posteriors[homGenotype.ordinal()] };
            allelePosteriors = MathUtils.normalizeFromLog10(allelePosteriors);

            // calculate the posterior weighted frequencies
            double[] HWvalues = getHardyWeinbergValues(f);
            double samplePofDgivenAFi = 0.0;
            samplePofDgivenAFi += HWvalues[GenotypeType.REF.ordinal()] * allelePosteriors[GenotypeType.REF.ordinal()];
            samplePofDgivenAFi += HWvalues[GenotypeType.HET.ordinal()] * allelePosteriors[GenotypeType.HET.ordinal()];
            samplePofDgivenAFi += HWvalues[GenotypeType.HOM.ordinal()] * allelePosteriors[GenotypeType.HOM.ordinal()];
            PofDgivenAFi += Math.log10(samplePofDgivenAFi);
        }

        return PofDgivenAFi;
    }

    protected List<Genotype> makeGenotypeCalls(char ref, HashMap<String, AlignmentContextBySample> contexts, GenomeLoc loc) {
        ArrayList<Genotype> calls = new ArrayList<Genotype>();

        for ( String sample : GLs.keySet() ) {

            // create the call
            Genotype call = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT, ref, loc);

            if ( call instanceof ReadBacked ) {
                ReadBackedPileup pileup = new ReadBackedPileup(ref, contexts.get(sample).getContext(StratifiedContext.OVERALL));
                ((ReadBacked)call).setPileup(pileup);
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(GLs.get(sample).getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(GLs.get(sample).getPosteriors());
            }

            calls.add(call);
        }

        return calls;
    }
}