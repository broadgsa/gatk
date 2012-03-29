package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager.SetOfStates;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

public abstract class VariantStratifier implements Comparable<VariantStratifier>, SetOfStates {
    private VariantEvalWalker variantEvalWalker;
    final private String name;
    final protected ArrayList<Object> states = new ArrayList<Object>();

    protected VariantStratifier() {
        name = this.getClass().getSimpleName();
    }

    // -------------------------------------------------------------------------------------
    //
    // to be overloaded
    //
    // -------------------------------------------------------------------------------------

    public abstract void initialize();

    public abstract List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName);

    // -------------------------------------------------------------------------------------
    //
    // final capabilities
    //
    // -------------------------------------------------------------------------------------

    /**
     * @return a reference to the parent VariantEvalWalker running this stratification
     */
    public final VariantEvalWalker getVariantEvalWalker() {
        return variantEvalWalker;
    }

    /**
     * Should only be called by VariantEvalWalker itself
     * @param variantEvalWalker
     */
    public final void setVariantEvalWalker(VariantEvalWalker variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
    }

    public final int compareTo(VariantStratifier o1) {
        return this.getName().compareTo(o1.getName());
    }

    public final String getName() {
        return name;
    }

    public final ArrayList<Object> getAllStates() {
        return states;
    }
}
