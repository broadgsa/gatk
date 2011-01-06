package org.broadinstitute.sting.playground.gatk.walkers.newvarianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;

import java.util.ArrayList;
import java.util.Set;

public class FunctionalClassStratifier extends VariantStratifier {
    // needs to know the variant context
    private ArrayList<String> states;

    @Override
    public void initialize(Set<VariantContextUtils.JexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames) {
        states = new ArrayList<String>();
        states.add("all");
        states.add("silent");
        states.add("missense");
        states.add("nonsense");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, VariantContext comp, VariantContext eval, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            String type = null;

            if (eval.getAttributeAsString("refseq.functionalClass") != null) {
                type = eval.getAttributeAsString("refseq.functionalClass");
            } else if (eval.getAttributeAsString("refseq.functionalClass_1") != null) {
                int annotationId = 1;
                String key;

                do {
                    key = String.format("refseq.functionalClass_%d", annotationId);

                    String newtype = eval.getAttributeAsString(key);

                    if ( newtype != null &&
                         ( type == null ||
                         ( type.equals("silent") && !newtype.equals("silent") ) ||
                         ( type.equals("missense") && newtype.equals("nonsense") ) )
                       ) {
                        type = newtype;
                    }

                    annotationId++;
                } while (eval.getAttributeAsString(key) != null);
            }

            if (type != null) {
                if      (type.equals("silent"))   { relevantStates.add("silent");   }
                else if (type.equals("missense")) { relevantStates.add("missense"); }
                else if (type.equals("nonsense")) { relevantStates.add("nonsense"); }
            }
        }

        return relevantStates;
    }
}
