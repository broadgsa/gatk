package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies the evaluation by each contig in the reference sequence
 */
public class Contig extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getContigNames());
        states.add("all");
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        if (eval != null) {
            relevantStates.add("all");
            relevantStates.add(eval.getChr());
        }

        return relevantStates;
    }
}
