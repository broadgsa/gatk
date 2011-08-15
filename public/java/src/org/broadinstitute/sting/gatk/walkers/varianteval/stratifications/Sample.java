package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;

public class Sample extends VariantStratifier {
    // needs the sample names
    private ArrayList<String> samples;

    @Override
    public void initialize() {
        samples = new ArrayList<String>();
        samples.addAll(getVariantEvalWalker().getSampleNamesForEvaluation());
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
