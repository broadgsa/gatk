package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;

import java.util.Set;


public class SingleSampleGenotypeCalculationModel extends GenotypeCalculationModel {

    private CallResult callsMetrics;

    public SingleSampleGenotypeCalculationModel(BaseMismatchModel baseModel,
                                                Set<String> samples,
                                                EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform platform,
                                                GenotypeWriter out,
                                                boolean genotypeMode,
                                                double lod,
                                                int maxDeletions,
                                                boolean verbose) {
        super(baseModel, samples, platform, out, genotypeMode, lod, maxDeletions, verbose);

        // check that this truly is a single-sample bam
        if ( samples.size() > 1 )
            throw new StingException("Single-sample genotype calculation model used on multi-sample bam... aborting");

        callsMetrics = new CallResult();
    }

    /**
     * Cloning of the object
     * @return clone
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public boolean calculateGenotype(char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        GenotypeLikelihoods gl = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
        gl.setVerbose(VERBOSE);
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);

        // check that there aren't too many deletions in the pileup
        pileup.setIncludeDeletionsInPileupString(true);
        if ( Utils.countOccurrences(BasicPileup.DELETION_CHAR, pileup.getBases()) > maxDeletionsInPileup )
            return false;

        gl.add(pileup, true);
        gl.validate();
        SSGenotypeCall call = new SSGenotypeCall(context.getLocation(), ref, gl, pileup);

        callsMetrics.nCalledBases++;

        if ( GENOTYPE_MODE || call.isVariant(call.getReference()) ) {
            double confidence = (GENOTYPE_MODE ? call.getNegLog10PError() : call.toVariation().getNegLog10PError());
            if ( confidence >= LOD_THRESHOLD ) {
                callsMetrics.nConfidentCalls++;
                //System.out.printf("Call %s%n", call);
                out.addGenotypeCall(call);
            } else
                callsMetrics.nNonConfidentCalls++;
        }
        return true;
    }

    public class CallResult {
        long nConfidentCalls = 0;
        long nNonConfidentCalls = 0;
        long nCalledBases = 0;

        CallResult() {}

        public String toString() {
            return String.format("SSG: %d confident and %d non-confident calls were made at %d bases",
                    nConfidentCalls, nNonConfidentCalls, nCalledBases);
        }
    }
}
