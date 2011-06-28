package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;

import java.util.ArrayList;
import java.util.Set;

public class Contig extends VariantStratifier {
    // needs to know the variant context
    private ArrayList<String> states;

    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames) {
        states = new ArrayList<String>();
        states.addAll(contigNames);
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
