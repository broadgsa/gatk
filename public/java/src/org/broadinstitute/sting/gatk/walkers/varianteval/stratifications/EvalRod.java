package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

public class EvalRod extends VariantStratifier implements RequiredStratification {
    private ArrayList<String> states;

    @Override
    public void initialize() {
        states = new ArrayList<String>();
        for ( RodBinding<VariantContext> rod : getVariantEvalWalker().getEvals() )
            states.add(rod.getName());
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add(evalName);

        return relevantStates;
    }
}
