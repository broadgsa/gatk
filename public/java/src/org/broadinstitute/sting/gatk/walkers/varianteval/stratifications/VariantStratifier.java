package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;

import java.util.ArrayList;
import java.util.Set;

public abstract class VariantStratifier implements Comparable {
    private VariantEvalWalker variantEvalWalker;

    /**
     * @return a reference to the parent VariantEvalWalker running this stratification
     */
    public VariantEvalWalker getVariantEvalWalker() {
        return variantEvalWalker;
    }

    /**
     * Should only be called by VariantEvalWalker itself
     * @param variantEvalWalker
     */
    public void setVariantEvalWalker(VariantEvalWalker variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
    }

    public abstract void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames);

    public ArrayList<String> getAllStates() {
        return new ArrayList<String>();
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        return null;
    }

    public int compareTo(Object o1) {
        return this.getClass().getSimpleName().compareTo(o1.getClass().getSimpleName());
    }
}
