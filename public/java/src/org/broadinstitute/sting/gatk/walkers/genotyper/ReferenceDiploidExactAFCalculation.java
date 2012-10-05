package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;

public class ReferenceDiploidExactAFCalculation extends DiploidExactAFCalculation {
    public ReferenceDiploidExactAFCalculation(final int nSamples, final int maxAltAlleles) {
        super(nSamples, maxAltAlleles);
    }

    public ReferenceDiploidExactAFCalculation(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
    }

    protected MaxLikelihoodSeen makeMaxLikelihood(final VariantContext vc, final AlleleFrequencyCalculationResult result) {
        return new ExactAFCalculation.MaxLikelihoodSeen();
    }
}
