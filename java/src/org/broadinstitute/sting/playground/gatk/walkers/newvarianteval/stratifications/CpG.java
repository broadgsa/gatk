package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Set;

public class CpG extends VariantStratifier {
    private ArrayList<String> states;

    @Override
    public void initialize(Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames) {
        states = new ArrayList<String>();
        states.add("all");
        states.add("CpG");
        states.add("non_CpG");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, String compName, VariantContext eval, String sampleName) {
        boolean isCpG = false;
        if (ref != null && ref.getBases() != null) {
            String fwRefBases = new String(ref.getBases());

            String leftFlank = fwRefBases.substring((fwRefBases.length()/2) - 1, (fwRefBases.length()/2) + 1);
            String rightFlank = fwRefBases.substring((fwRefBases.length()/2), (fwRefBases.length()/2) + 2);

            if (leftFlank.equalsIgnoreCase("CG") || leftFlank.equalsIgnoreCase("GC") || rightFlank.equalsIgnoreCase("CG") || rightFlank.equalsIgnoreCase("GC")) {
                isCpG = true;
            }
        }

        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add("all");
        relevantStates.add(isCpG ? "CpG" : "non_CpG");

        return relevantStates;
    }
}
