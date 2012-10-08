package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;

public class ReferenceDiploidExactAFCalc extends DiploidExactAFCalc {
    public ReferenceDiploidExactAFCalc(final int nSamples, final int maxAltAlleles) {
        super(nSamples, maxAltAlleles);
    }

    public ReferenceDiploidExactAFCalc(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
    }

    protected StateTracker makeMaxLikelihood(final VariantContext vc, final AFCalcResultTracker resultTracker) {
        return new StateTracker();
    }
}
