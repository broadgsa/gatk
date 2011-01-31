package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Set;

public class JexlExpression extends VariantStratifier implements StandardStratification {
    // needs to know the jexl expressions
    private Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions;
    private ArrayList<String> states;

    @Override
    public void initialize(Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames) {
        this.jexlExpressions = jexlExpressions;

        states = new ArrayList<String>();
        states.add("none");
        for ( VariantContextUtils.JexlVCMatchExp jexlExpression : jexlExpressions ) {
            states.add(jexlExpression.name);
        }
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, String compName, VariantContext eval, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add("none");

        for ( VariantContextUtils.JexlVCMatchExp jexlExpression : jexlExpressions ) {
            if (eval != null && VariantContextUtils.match(eval, jexlExpression)) {
                relevantStates.add(jexlExpression.name);
            }
        }

        return relevantStates;
    }
}
