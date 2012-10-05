package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;

public class ConstrainedDiploidExactAFCalculation extends DiploidExactAFCalculation {
    public ConstrainedDiploidExactAFCalculation(final int nSamples, final int maxAltAlleles) {
        super(nSamples, maxAltAlleles);
    }

    public ConstrainedDiploidExactAFCalculation(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
    }

    protected MaxLikelihoodSeen makeMaxLikelihood(final VariantContext vc, final AlleleFrequencyCalculationResult result) {
        final int[] maxACsToConsider = computeMaxACs(vc);
        result.setAClimits(maxACsToConsider);
        return new MaxLikelihoodSeen(maxACsToConsider);
    }
}
