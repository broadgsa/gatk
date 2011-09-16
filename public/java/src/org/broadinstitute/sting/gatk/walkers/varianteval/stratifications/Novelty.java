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


    @Override
    public void initialize() {
        states = new ArrayList<String>(Arrays.asList("all", "known", "novel"));
        knowns = getVariantEvalWalker().getKnowns();
    }

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (tracker != null && eval != null) {
            final Collection<VariantContext> knownComps = tracker.getValues(knowns, ref.getLocus());
            for ( final VariantContext c : knownComps ) {
                // loop over sites, looking for something that matches the type eval
                if ( eval.getType() == c.getType() ) {
                    return Arrays.asList("all", "known");
                }
            }
        }

        return Arrays.asList("all", "novel");
    }
}
