/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager.Stratifier;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public abstract class VariantStratifier implements Comparable<VariantStratifier>, Stratifier {
    private VariantEval variantEvalWalker;
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

    public abstract List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String familyName);

    // -------------------------------------------------------------------------------------
    //
    // final capabilities
    //
    // -------------------------------------------------------------------------------------

    /**
     * @return a reference to the parent VariantEvalWalker running this stratification
     */
    public final VariantEval getVariantEvalWalker() {
        return variantEvalWalker;
    }

    /**
     * Should only be called by VariantEvalWalker itself
     * @param variantEvalWalker
     */
    public final void setVariantEvalWalker(VariantEval variantEvalWalker) {
        this.variantEvalWalker = variantEvalWalker;
    }

    public final int compareTo(VariantStratifier o1) {
        return this.getName().compareTo(o1.getName());
    }

    @Override
    public String toString() {
        return getName();
    }

    public final String getName() {
        return name;
    }
    
    public String getFormat() { return "%s"; }
    
    public final ArrayList<Object> getAllStates() {
        return states;
    }


    /**
     * The way for a stratifier to specify that it's incompatible with specific evaluations.  For
     * example, VariantSummary includes a per-sample metric, and so cannot be used safely with Sample
     * or AlleleCount stratifications as this introduces an O(n^2) memory and cpu cost.
     *
     * @return the set of VariantEvaluators that cannot be active with this Stratification
     */
    public Set<Class<? extends VariantEvaluator>> getIncompatibleEvaluators() {
        return Collections.emptySet();
    }
}
