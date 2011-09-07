package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies the eval RODs by the allele frequency of the alternate allele
 *
 * Uses a constant 0.005 frequency grid, and projects the AF INFO field value.  Requires
 * that AF be present in every ROD, otherwise this stratification throws an exception
 */
public class AlleleFrequency extends VariantStratifier {
    @Override
    public void initialize() {
        states = new ArrayList<String>();
        for( double a = 0.000; a <= 1.005; a += 0.005 ) {
            states.add(String.format("%.3f", a));
        }
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        if (eval != null) {
            try {
                relevantStates.add(String.format("%.3f", (5.0 * MathUtils.round(eval.getAttributeAsDouble("AF", 0.0) / 5.0, 3))));
            } catch (Exception e) {
                return relevantStates;
            }
        }

        return relevantStates;
    }
}
