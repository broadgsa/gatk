package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public abstract class EMGenotypeCalculationModel extends GenotypeCalculationModel {

    // We need to set a limit on the EM iterations in case something flukey goes on
    protected static final int MAX_EM_ITERATIONS = 8;

    // We consider the EM stable when the MAF doesn't change more than 1/10,000
    protected static final double EM_STABILITY_METRIC = 1e-4;

    protected EMGenotypeCalculationModel() {}

    public VariantCallContext callLocus(RefMetaDataTracker tracker, char ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {

        // run the EM calculation
        EMOutput overall = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

        double PofD = Math.pow(10, overall.getPofD());
        double PofNull = Math.pow(10, overall.getPofNull());
        double sum = PofD + PofNull;

        // calculate the phred-scaled confidence score
        double phredScaledConfidence;
        boolean bestIsRef = false;
        if ( PofD > PofNull ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofNull / sum);
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofD / sum);
            bestIsRef = true;
        }

        // are we above the lod threshold for emitting calls (and not in all-bases mode)?
        if ( !ALL_BASE_MODE && ((!GENOTYPE_MODE && bestIsRef) || phredScaledConfidence < CONFIDENCE_THRESHOLD) )
            return new VariantCallContext(phredScaledConfidence >= CONFIDENCE_THRESHOLD);

        // generate the calls
        List<Genotype> calls = genotypeCallsFromGenotypeLikelihoods(overall, ref, contexts);

        VariationCall locusdata = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT, ref, loc, bestIsRef ? Variation.VARIANT_TYPE.REFERENCE : Variation.VARIANT_TYPE.SNP);
        if ( locusdata != null ) {
            if ( locusdata instanceof ConfidenceBacked ) {
                ((ConfidenceBacked)locusdata).setConfidence(phredScaledConfidence);
            }
            if ( locusdata instanceof IDBacked ) {
                rodDbSNP dbsnp = getDbSNP(tracker);
                if ( dbsnp != null )
                    ((IDBacked)locusdata).setID(dbsnp.getRS_ID());
            }
            if ( locusdata instanceof SLODBacked && REPORT_SLOD ) {

                // calculate strand score

                double lod = overall.getPofD() - overall.getPofNull();
                logger.debug(String.format("LOD=%f", lod));

                EMOutput forward = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.FORWARD);
                EMOutput reverse = runEM(ref, contexts, priors, StratifiedAlignmentContext.StratifiedContextType.REVERSE);
                double forwardLod = (forward.getPofD() + reverse.getPofNull()) - overall.getPofNull();
                double reverseLod = (reverse.getPofD() + forward.getPofNull()) - overall.getPofNull();
                logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);
                double strandScore = Math.max(forwardLod - lod, reverseLod - lod);
                logger.debug(String.format("SLOD=%f", strandScore));
                // rescale by a factor of 10
                strandScore *= 10.0;

                ((SLODBacked)locusdata).setSLOD(strandScore);
            }
            locusdata.setNonRefAlleleFrequency(overall.getMAF());

            // finally, associate the Variation with the Genotypes
            locusdata.setGenotypeCalls(calls);
            for ( Genotype call : calls )
                ((GenotypeCall)call).setVariation(locusdata);
        }
        return new VariantCallContext(locusdata, calls, phredScaledConfidence >= CONFIDENCE_THRESHOLD);
    }

    protected List<Genotype> genotypeCallsFromGenotypeLikelihoods(EMOutput results, char ref, Map<String, StratifiedAlignmentContext> contexts) {
        HashMap<String, GenotypeLikelihoods> GLs = results.getGenotypeLikelihoods();

        ArrayList<Genotype> calls = new ArrayList<Genotype>();
        int variantCalls = 0;

        for ( String sample : GLs.keySet() ) {

            // create the call
            AlignmentContext context = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);
            GenotypeCall call = GenotypeWriterFactory.createSupportedGenotypeCall(OUTPUT_FORMAT, ref, context.getLocation());

            // set the genotype and confidence
            double[] posteriors = GLs.get(sample).getPosteriors();
            Integer sorted[] = Utils.SortPermutation(posteriors);
            DiploidGenotype bestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];
            DiploidGenotype nextGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 2]];
            call.setNegLog10PError(posteriors[bestGenotype.ordinal()] - posteriors[nextGenotype.ordinal()]);
            call.setGenotype(bestGenotype);

            if ( call instanceof ReadBacked ) {
                ReadBackedPileup pileup = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
                ((ReadBacked)call).setPileup(pileup);
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(GLs.get(sample).getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(posteriors);
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

    public EMOutput runEM(char ref, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType) {

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
    protected abstract void initializeGenotypeLikelihoods(char ref, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors, StratifiedAlignmentContext.StratifiedContextType contextType);
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