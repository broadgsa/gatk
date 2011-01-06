package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Set;

public class JexlExpressionStratifier extends VariantStratifier implements StandardStratification {
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

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, VariantContext eval, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add("none");

        for ( VariantContextUtils.JexlVCMatchExp jexlExpression : jexlExpressions ) {
            System.out.println(jexlExpression.name + " " + jexlExpression.exp.getExpression());
            if (VariantContextUtils.match(eval, jexlExpression)) {
                relevantStates.add(jexlExpression.name);
            }
        }

        return relevantStates;
    }
}
