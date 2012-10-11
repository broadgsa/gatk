package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public class ReferenceDiploidExactAFCalc extends DiploidExactAFCalc {
    protected ReferenceDiploidExactAFCalc(int nSamples, int maxAltAlleles, int maxAltAllelesForIndels, final int ploidy) {
        super(nSamples, maxAltAlleles, maxAltAllelesForIndels, ploidy);
    }

    protected StateTracker makeMaxLikelihood(final VariantContext vc, final AFCalcResultTracker resultTracker) {
        return new StateTracker();
    }
}
