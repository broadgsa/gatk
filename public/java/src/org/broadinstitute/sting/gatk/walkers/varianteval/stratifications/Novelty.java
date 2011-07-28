package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.Set;

public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private Set<String> knownNames;
    private ArrayList<String> states;

    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames) {
        this.knownNames = knownNames;

        states = new ArrayList<String>();
        states.add("all");
        states.add("known");
        states.add("novel");
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        boolean isNovel = true;

        if (tracker != null) {
            for (String knownName : knownNames) {
                if (tracker.hasValues(knownName)) {
                    EnumSet<VariantContext.Type> allowableTypes = EnumSet.of(VariantContext.Type.NO_VARIATION);
                    if (eval != null) {
                        allowableTypes.add(eval.getType());
                    }

                    Collection<VariantContext> knownComps = tracker.getVariantContexts(ref, knownName, allowableTypes, ref.getLocus(), true, true);

                    isNovel = knownComps.size() == 0;

                    break;
                }
            }
        }

        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add("all");
        relevantStates.add(isNovel ? "novel" : "known");

        return relevantStates;
    }
}
