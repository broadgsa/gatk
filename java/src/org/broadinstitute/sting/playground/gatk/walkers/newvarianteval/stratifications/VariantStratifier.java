package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.ArrayList;
import java.util.Set;

public abstract class VariantStratifier implements Comparable {
    public abstract void initialize(Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames);

    public ArrayList<String> getAllStates() {
        return new ArrayList<String>();
    }

    public boolean isApplicable(RefMetaDataTracker tracker, ReferenceContext ref, VariantContext comp, VariantContext eval, String state) {
        return false;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, VariantContext eval, String sampleName) {
        return null;
    }

    public int compareTo(Object o1) {
        return this.getClass().getSimpleName().compareTo(o1.getClass().getSimpleName());
    }
}
