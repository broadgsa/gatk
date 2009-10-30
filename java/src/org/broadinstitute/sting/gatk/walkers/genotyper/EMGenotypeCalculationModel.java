package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public abstract class EMGenotypeCalculationModel extends GenotypeCalculationModel {

    // We need to set a limit on the EM iterations in case something flukey goes on
    protected static final int MAX_EM_ITERATIONS = 8;

    // We consider the EM stable when the MAF doesn't change more than 1/10,000
    protected static final double EM_STABILITY_METRIC = 1e-4;

    protected EMGenotypeCalculationModel() {}

    public Pair<List<Genotype>, GenotypeMetaData> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the context for each sample, overall and separated by strand
        HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context);
        if ( contexts == null )
            return null;

        // run the EM calculation
        EMOutput overall = runEM(ref, contexts, priors, StratifiedContext.OVERALL);
        double lod = overall.getPofD() - overall.getPofNull();

        // return a null call if we don't pass the lod cutoff
        if ( !ALL_BASE_MODE && lod < LOD_THRESHOLD )
            return new Pair<List<Genotype>, GenotypeMetaData>(null, null);

        // calculate strand score
        EMOutput forward = runEM(ref, contexts, priors, StratifiedContext.FORWARD);
        EMOutput reverse = runEM(ref, contexts, priors, StratifiedContext.REVERSE);
        double forwardLod = (forward.getPofD() + reverse.getPofNull()) - overall.getPofNull();
        double reverseLod = (reverse.getPofD() + forward.getPofNull()) - overall.getPofNull();
        logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);
        double strandScore = Math.max(forwardLod - lod, reverseLod - lod);

        logger.debug(String.format("LOD=%f, SLOD=%f", lod, strandScore));

        // generate the calls
        GenotypeMetaData metadata = new GenotypeMetaData(lod, strandScore, overall.getMAF());
        return new Pair<List<Genotype>, GenotypeMetaData>(genotypeCallsFromGenotypeLikelihoods(overall, ref, contexts), metadata);
    }

    protected List<Genotype> genotypeCallsFromGenotypeLikelihoods(EMOutput results, char ref, HashMap<String, AlignmentContextBySample> contexts) {
        HashMap<String, GenotypeLikelihoods> GLs = results.getGenotypeLikelihoods();

        ArrayList<Genotype> calls = new ArrayList<Genotype>();
        int variantCalls = 0;

        for ( String sample : GLs.keySet() ) {

            // create the call
            AlignmentContext context = contexts.get(sample).getContext(StratifiedContext.OVERALL);
            Genotype call = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT, ref, context.getLocation());

            if ( call instanceof ReadBacked ) {
                ((ReadBacked)call).setReads(context.getReads());
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

            // increment the variant count if it's non-ref
            if ( call.isVariant(ref) )
                variantCalls++;
        }

        // if everyone is ref, don't emit any calls
        if ( variantCalls == 0 )
            calls.clear();

        return calls;
    }

    public EMOutput runEM(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {

        // initialize the allele frequencies
        initializeAlleleFrequencies(contexts.size(), ref);

        // initialize the genotype likelihoods
        initializeGenotypeLikelihoods(ref, contexts, priors, contextType);

        // The EM loop:
        //   we want to continue until the calculation is stable, but we need some max on the number of iterations
        int iterations = 0;
        boolean EM_IS_STABLE;

        do {
            calculateAlleleFrequencyPosteriors();

            applyAlleleFrequencyToGenotypeLikelihoods();

            EM_IS_STABLE = isStable();

        } while ( ++iterations < MAX_EM_ITERATIONS && !EM_IS_STABLE );

        logger.debug("EM loop took " + iterations + " iterations");

        return computePofF(ref);
    }

    protected abstract void initializeAlleleFrequencies(int numSamplesInContext, char ref);
    protected abstract void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType);
    protected abstract void calculateAlleleFrequencyPosteriors();
    protected abstract void applyAlleleFrequencyToGenotypeLikelihoods();
    protected abstract boolean isStable();
    protected abstract EMOutput computePofF(char ref);


    /**
     * A class to keep track of the EM output
     */
    protected class EMOutput {
        private double pD, pNull, pF, MAF;
        private HashMap<String, GenotypeLikelihoods> GLs;

        EMOutput(double pD, double pNull, double pF, double MAF, HashMap<String, GenotypeLikelihoods> GLs) {
            this.pD = pD;
            this.pNull = pNull;
            this.pF = pF;
            this.MAF = MAF;
            this.GLs = GLs;
        }

        public double getPofD() { return pD; }
        public double getPofNull() { return pNull; }
        public double getPofF() { return pF; }
        public double getMAF() { return MAF; }
        public HashMap<String, GenotypeLikelihoods> getGenotypeLikelihoods() { return GLs; }
    }
}