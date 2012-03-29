package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Stratifies the evaluation by each contig in the reference sequence
 */
public class Contig extends VariantStratifier {
    @Override
    public void initialize() {
        states.addAll(getVariantEvalWalker().getContigNames());
        states.add("all");
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (eval != null) {
            return Arrays.asList((Object)"all", eval.getChr());
        } else {
            return Collections.emptyList();
        }
    }
}
