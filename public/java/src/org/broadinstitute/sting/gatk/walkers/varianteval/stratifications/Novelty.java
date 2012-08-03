package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * Stratifies by whether a site in in the list of known RODs (e.g., dbsnp by default)
 */
public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private List<RodBinding<VariantContext>> knowns;

    private final static List<Object> KNOWN_STATES = Arrays.asList((Object)"all", (Object)"known");
    private final static List<Object> NOVEL_STATES = Arrays.asList((Object)"all", (Object)"novel");

    @Override
    public void initialize() {
        states.addAll(Arrays.asList("all", "known", "novel"));
        knowns = getVariantEvalWalker().getKnowns();
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (tracker != null && eval != null) {
            final Collection<VariantContext> knownComps = tracker.getValues(knowns, ref.getLocus());
            for ( final VariantContext c : knownComps ) {
                // loop over sites, looking for something that matches the type eval
                if ( eval.getType() == c.getType() || eval.getType() == VariantContext.Type.NO_VARIATION ) {
                    return KNOWN_STATES;
                }
            }
        } 
        
        return NOVEL_STATES;
    }
}