package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

public class Filter extends VariantStratifier {
    // needs to know the variant context
    private ArrayList<String> states;

    @Override
    public void initialize() {
        states = new ArrayList<String>();
        states.add("called");
        states.add("filtered");
        states.add("raw");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("raw");
        if (eval != null) {
            relevantStates.add(eval.isFiltered() ? "filtered" : "called");
        }

        return relevantStates;
    }
}
