package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public abstract class VariantStratifier implements Comparable<VariantStratifier> {
    private VariantEvalWalker variantEvalWalker;
    final private String name;
    protected ArrayList<String> states = new ArrayList<String>();

    protected VariantStratifier() {
        name = this.getClass().getSimpleName();
    }

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

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        return null;
    }

    public int compareTo(VariantStratifier o1) {
        return this.getName().compareTo(o1.getName());
    }

    public final String getName() {
        return name;
    }

    public ArrayList<String> getAllStates() {
        return states;
    }
}
