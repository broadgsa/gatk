package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Set;

public class Sample extends VariantStratifier {
    // needs the sample names
    private ArrayList<String> samples;

    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames) {
        samples = new ArrayList<String>();
        samples.addAll(sampleNames);
    }

    public ArrayList<String> getAllStates() {
        return samples;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add(sampleName);

        return relevantStates;
    }
}
