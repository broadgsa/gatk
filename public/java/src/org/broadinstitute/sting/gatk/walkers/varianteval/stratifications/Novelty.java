package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private Set<String> knownNames;
    final private ArrayList<String> states = new ArrayList<String>(Arrays.asList("all", "known", "novel"));


    @Override
    public void initialize(Set<SortableJexlVCMatchExp> jexlExpressions, Set<String> compNames, Set<String> knownNames, Set<String> evalNames, Set<String> sampleNames, Set<String> contigNames) {
        this.knownNames = knownNames;
    }

    public ArrayList<String> getAllStates() {
        return states;
    }

    public ArrayList<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (tracker != null && eval != null) {
            for (final String knownName : knownNames) {
                final Collection<VariantContext> knownComps = tracker.getVariantContexts(knownName, ref.getLocus(), true, false);
                for ( final VariantContext c : knownComps ) {
                    // loop over sites, looking for something that matches the type eval
                    if ( eval.getType() == c.getType() ) {
                        return new ArrayList<String>(Arrays.asList("all", "known"));
                    }
                }
            }
        }

        return new ArrayList<String>(Arrays.asList("all", "novel"));
    }
}
