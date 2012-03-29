package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collections;
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
        for( double a = 0.000; a <= 1.005; a += 0.005 ) {
            states.add(String.format("%.3f", a));
        }
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (eval != null) {
            try {
                return Collections.singletonList((Object)String.format("%.3f", (5.0 * MathUtils.round(eval.getAttributeAsDouble("AF", 0.0) / 5.0, 3))));
            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }
}
