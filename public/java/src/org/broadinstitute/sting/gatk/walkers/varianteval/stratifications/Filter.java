package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by the FILTER status (PASS, FAIL) of the eval records
 */
public class Filter extends VariantStratifier {
    @Override
    public void initialize() {
        states.add("called");
        states.add("filtered");
        states.add("raw");
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("raw");
        if (eval != null) {
            relevantStates.add(eval.isFiltered() ? "filtered" : "called");
        }

        return relevantStates;
    }
}
