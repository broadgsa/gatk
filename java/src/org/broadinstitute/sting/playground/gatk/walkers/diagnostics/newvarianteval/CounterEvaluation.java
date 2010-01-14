package org.broadinstitute.sting.playground.gatk.walkers.diagnostics.newvarianteval;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.List;
import java.util.ArrayList;

public class CounterEvaluation implements VariantEvaluation {
    private int numVariants = 0;
    private int numCalls = 0;
    private int numKnown = 0;

    public CounterEvaluation() {}

    public void update(RodVCF variant, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        numVariants += (variant != null && variant.isSNP()) ? 1 : 0;
        numKnown += (variant != null && variant.isSNP() && !variant.getID().equals(".")) ? 1 : 0;
        numCalls++;
    }

    public String getName() {
        return "Counter";
    }

    public List<String> getResult() {
        ArrayList<String> results = new ArrayList<String>();

        double variantRate = ((double) numVariants)/((double) numCalls);
        double knownVariantsPct = 100.0*((double) numKnown)/((double) numVariants);

        results.add(String.format("regionSize=%d", numCalls));
        results.add(String.format("variants=%d", numVariants));
        results.add(String.format("variantsKnown=%d", numKnown));
        results.add(String.format("knownVariantsPct=%f", knownVariantsPct));
        results.add(String.format("variantRate=1/%d", Math.round(1.0/variantRate)));

        return results;
    }
}
