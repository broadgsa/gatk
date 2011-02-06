package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class Degeneracy extends VariantStratifier {
    private ArrayList<String> states;

    private HashMap<String, String> degeneracies;

    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames) {
        states = new ArrayList<String>();
        states.add("1-fold");
        states.add("2-fold");
        states.add("3-fold");
        states.add("4-fold");
        states.add("6-fold");
        states.add("all");

        degeneracies = new HashMap<String, String>();
        degeneracies.put("Ile", "3-fold");
        degeneracies.put("Leu", "6-fold");
        degeneracies.put("Val", "4-fold");
        degeneracies.put("Phe", "2-fold");
        degeneracies.put("Met", "1-fold");
        degeneracies.put("Cys", "2-fold");
        degeneracies.put("Ala", "4-fold");
        degeneracies.put("Gly", "4-fold");
        degeneracies.put("Pro", "4-fold");
        degeneracies.put("Thr", "4-fold");
        degeneracies.put("Ser", "6-fold");
        degeneracies.put("Tyr", "2-fold");
        degeneracies.put("Try", "1-fold");
        degeneracies.put("Trp", "1-fold");
        degeneracies.put("Gln", "2-fold");
        degeneracies.put("Asn", "2-fold");
        degeneracies.put("His", "2-fold");
        degeneracies.put("Glu", "2-fold");
        degeneracies.put("Asp", "2-fold");
        degeneracies.put("Lys", "2-fold");
        degeneracies.put("Arg", "6-fold");
        degeneracies.put("Stop", "3-fold");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            String type = null;
            String aa = null;

            if (eval.hasAttribute("refseq.functionalClass")) {
                aa = eval.getAttributeAsString("refseq.variantAA");
            } else if (eval.hasAttribute("refseq.functionalClass_1")) {
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

                        String aakey = String.format("refseq.variantAA_%d", annotationId);
                        aa = eval.getAttributeAsString(aakey);
                    }

                    annotationId++;
                } while (eval.hasAttribute(key));
            }

            if (aa != null && degeneracies.containsKey(aa)) {
                relevantStates.add(degeneracies.get(aa));
            }
        }

        return relevantStates;
    }
}
