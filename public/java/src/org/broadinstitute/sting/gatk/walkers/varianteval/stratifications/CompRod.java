package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;


/**
 * Required stratification grouping output by each comp ROD
 */
public class CompRod extends VariantStratifier implements RequiredStratification {
    @Override
    public void initialize() {
        for ( RodBinding<VariantContext> rod : getVariantEvalWalker().getComps() )
            states.add(rod.getName());
    }


    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add(compName);

        return relevantStates;
    }
}
