package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the eval RODs by each sample in the eval ROD.
 *
 * This allows the system to analyze each sample separately.  Since many evaluations
 * only consider non-reference sites, stratifying by sample results in meaningful
 * calculations for CompOverlap
 */
public class Sample extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getSampleNamesForStratification());
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        return Arrays.asList(sampleName);
    }
}
