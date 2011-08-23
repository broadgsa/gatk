package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

public class Contig extends VariantStratifier {
    // needs to know the variant context
    private ArrayList<String> states;

    @Override
    public void initialize() {
        states = new ArrayList<String>();
        states.addAll(getVariantEvalWalker().getContigNames());
        states.add("all");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        if (eval != null) {
            relevantStates.add("all");
            relevantStates.add(eval.getChr());
        }

        return relevantStates;
    }
}
