package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;

import java.util.ArrayList;
import java.util.Set;

public class FunctionalClass extends VariantStratifier {
    // needs to know the variant context
    private ArrayList<String> states;

    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames) {
        states = new ArrayList<String>();
        states.add("all");
        states.add("silent");
        states.add("missense");
        states.add("nonsense");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            String type = null;

            if (eval.hasAttribute("refseq.functionalClass")) {
                type = eval.getAttributeAsString("refseq.functionalClass");
            } else if (eval.hasAttribute("refseq.functionalClass_1")) {
                int annotationId = 1;
                String key;

                do {
                    key = String.format("refseq.functionalClass_%d", annotationId);

                    String newtype = eval.getAttributeAsString(key);

                    if ( newtype != null && !newtype.equalsIgnoreCase("null") &&
                         ( type == null ||
                         ( type.equals("silent") && !newtype.equals("silent") ) ||
                         ( type.equals("missense") && newtype.equals("nonsense") ) )
                       ) {
                        type = newtype;
                    }

                    annotationId++;
                } while (eval.hasAttribute(key));
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
