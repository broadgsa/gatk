package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Set;

public class CompRod extends VariantStratifier implements RequiredStratification {
    // Needs to know the comp rods
    private Set<String> compNames;
    private ArrayList<String> states;

    @Override
    public void initialize(Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames) {
        this.compNames = compNames;

        states = new ArrayList<String>();
        states.addAll(compNames);
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, String compName, VariantContext eval, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add(compName);

        return relevantStates;
    }
}
