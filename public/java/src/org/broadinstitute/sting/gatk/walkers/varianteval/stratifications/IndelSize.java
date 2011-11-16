package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies the eval RODs by the indel size
 *
 * Indel sizes are stratified from sizes -100 to +100. Sizes greater than this are lumped in the +/- 100 bin
 * This stratification ignores multi-allelic indels (whose size is not defined uniquely)
 */
public class IndelSize extends VariantStratifier {
    static final int MAX_INDEL_SIZE = 100;
    @Override
    public void initialize() {
        states = new ArrayList<String>();
        for( int a=-MAX_INDEL_SIZE; a <=MAX_INDEL_SIZE; a++ ) {
            states.add(String.format("%d", a));
        }
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        if (eval != null && eval.isIndel() && eval.isBiallelic()) {
            try {
                int eventLength = 0;
                if ( eval.isSimpleInsertion() ) {
                    eventLength = eval.getAlternateAllele(0).length();
                } else if ( eval.isSimpleDeletion() ) {
                    eventLength = -eval.getReference().length();
                }

                if (eventLength > MAX_INDEL_SIZE)
                    eventLength = MAX_INDEL_SIZE;
                else if (eventLength < -MAX_INDEL_SIZE)
                    eventLength = -MAX_INDEL_SIZE;

                relevantStates.add(String.format("%d",eventLength));
            } catch (Exception e) {
                return relevantStates;
            }
        }

        return relevantStates;
    }
}
