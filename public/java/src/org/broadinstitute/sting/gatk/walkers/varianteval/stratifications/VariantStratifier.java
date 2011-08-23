package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

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

    public abstract void initialize();

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
